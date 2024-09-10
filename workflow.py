
from gwf import Workflow
from pathlib import Path
import os, string, re

#################################################################################
# Templates
#################################################################################

def vcf2bed(vcf_file, bed_file):

    options = {'memory': '4g',
               'walltime': '11:00:00' }
    spec = """

    # source activate hri
    python ./scripts/format_single_site_bed.py {vcf_file} {bed_file}

    """.format(vcf_file=vcf_file, bed_file=bed_file)

    return [str(vcf_file)], [str(bed_file)], options, spec


def stevison2bed(map_file, bed_file):

    options = {'memory': '4g',
               'walltime': '11:00:00' }

    spec = """

    # source activate hri
    python ./scripts/stevison2bed.py {map_file} {bed_file}

    """.format(map_file=map_file, bed_file=bed_file)

    return [str(map_file)], [str(bed_file)], options, spec


def intervals_to_sites(intervals_file, sites_file): 

    options = {'memory': '4g',
               'walltime': '11:00:00'
               }

    spec = """

    # source activate hri
    python ./scripts/liftover_sites.py {intervals_file} {sites_file}

    """

    return [str(intervals_file)], [str(sites_file)], options, spec


def liftover(bed_file, chain_file, mapped_file, unmapped_file):

    options = {'memory': '4g',
               'walltime': '36:00:00'
               }

    spec = """

    # source activate hri
    liftOver -bedPlus=3 {bed_file} {chain_file} {mapped_file} {unmapped_file}

    """.format(bed_file=bed_file, chain_file=chain_file, mapped_file=mapped_file, unmapped_file=unmapped_file)

    return [str(bed_file), str(chain_file)], [str(mapped_file), str(unmapped_file)], options, spec


def bed_difference(bed_file1, bed_file2, output_file): 

    options = {'memory': '4g',
               'walltime': '11:00:00'
               }
    shell_spec = """

    # source activate hri
    python ./scripts/liftover_diff.py {bed_file1} {bed_file2} | sort -k1,1 -k2,2n -k3,3n > {output_file}

    """.format(bed_file1=bed_file1, bed_file2=bed_file2, output_file=output_file)

    return [str(bed_file1), str(bed_file2)], [str(output_file)], options, shell_spec


def bed_merge_and_split(input_files, output_files):

    output_base_names = [Path(x).name for x in output_files]

    # get output dir and make sure it is unambiguous:
    output_dirs = [Path(x).parent for x in output_files]
    assert len(set(output_dirs)) == 1
    output_dir = output_dirs[0]

    input_files = [str(x) for x in input_files]
    output_files = [str(x) for x in output_files]

    options = {'memory': '10g',
               'walltime': '11:00:00'}

    shell_spec = """

    sort -k1,1 -k2,2n -k3,3n --merge {input_files} -T /scratch/$GWF_JOBID | python ./scripts/bed_split.py {output_dir} {output_base_names}

    """.format(input_files=" ".join(input_files),
               output_dir=output_dir,
               output_base_names=" ".join(output_base_names))

    return list(map(str, input_files)), list(map(str, output_files)), options, shell_spec


def interference_score(sites_file, annot_file, map_file, score_file, map_scale, rate_scale):

    options = {'memory': '10g',
               'walltime': '11:00:00' }

    shell_spec = """

    # source activate hri
    python ./scripts/interference_map.py --map-col-idx 0,1,2,4 --map-scale {map_scale} --rate-scale {rate_scale} \
           --output {score_file}  {map_file} {annot_file} {sites_file} 

    """.format(sites_file=sites_file,
               annot_file=annot_file,
               map_file=map_file,
               score_file=score_file,
               map_scale=map_scale,
               rate_scale=rate_scale)

    return [str(sites_file), str(annot_file), str(map_file)], [str(score_file)], options, shell_spec


#################################################################################
# Workflow components
#################################################################################

gwf = Workflow(defaults={'account': 'hri'})


def split_file(vcf_file, split_sites_files_dir, n_files=50):
    """
    Split a file into 20 chunks.
    """
    if not split_sites_files_dir.exists():
        os.makedirs(str(split_sites_files_dir))
    prefix = str(split_sites_files_dir / vcf_file.with_suffix('').name)
    sites_files = [Path('{}.{}{}'.format(prefix, x, vcf_file.suffix)) for x in range(n_files)]
    gwf.target('ic_split_sites_{}'.format(species), inputs=[str(vcf_file)], outputs=[str(p) for p in sites_files]) << """

        python ./scripts/file_chunks.py $((1 + `wc -l < {vcf_file}`/{n_files})) {prefix} {suffix} {vcf_file}

        """.format(vcf_file=str(vcf_file), prefix=prefix, suffix=vcf_file.suffix, n_files=n_files)
    return sites_files


def format_vcf_to_bed(vcf_files, species, steps_dir):
    """
    Formats a list of VCF files to BED bed files.
    """
    if not steps_dir.exists():
        os.makedirs(str(steps_dir))
    bed_files = [steps_dir / x.with_suffix('.bed').name for x in vcf_files]
    for i, (vcf_file, bed_file) in enumerate(zip(vcf_files, bed_files)):
        gwf.target_from_template('ic_vcf2bed_{}_{}'.format(species, i),
            vcf2bed(vcf_file=str(vcf_file), bed_file=str(bed_file)))

    return bed_files


def format_orig_map_files_to_bed(map_files, steps_dir):
    """
    Formats a list of Stevison et al map files to BED format.
    """
    if not steps_dir.exists():
        os.makedirs(str(steps_dir))
    #bed_files = [steps_dir / x.with_suffix('.bed').name for x in map_files]
    bed_files = list()
    for p in map_files:
        file_name = re.match(r'[^-]+', p.with_suffix('.bed').name).group(0).lower() + '.bed'
        bed_files.append(steps_dir / file_name)

    for i, (map_file, bed_file) in enumerate(zip(map_files, bed_files)):
        gwf.target_from_template('ic_stevison2bed_{}_{}'.format(species, i),
            stevison2bed(map_file=str(map_file), bed_file=str(bed_file)))
    return bed_files


def reciprocal_liftover(intervals_files, forwards_chain_file, backwards_chain_file, 
                        slurm_tag, steps_dir, target_chromosomes):
    """
    Does reciprocal lift over of a set of intervals.
    """

    if not steps_dir.exists():
        os.makedirs(str(steps_dir))

    # output files
    mapped_files= [steps_dir / x.with_suffix('.mapped').name for x in intervals_files]
    unmapped_files = [x.with_suffix('.unmapped') for x in mapped_files]
    backmapped_files = [x.with_suffix('.backmapped') for x in mapped_files]
    unbackmapped_files = [x.with_suffix('.nobackmapped') for x in mapped_files]
    filtered_files = [x.with_suffix('.filtered') for x in mapped_files]

    lifted_files = [steps_dir / 'sorted' / "{}.bed".format(x) for x in target_chromosomes]

    for i, chrom in enumerate(intervals_files):

        # lift over intervals
        gwf.target_from_template('{}_lift_{}'.format(slurm_tag, i),
            liftover(bed_file=intervals_files[i], chain_file=forwards_chain_file,
            mapped_file=mapped_files[i], unmapped_file=unmapped_files[i]))

        # lift back to orginal coordinates to ensure one to one correspondence
        gwf.target_from_template('{}_liftback_{}'.format(slurm_tag, i),
            liftover(bed_file=mapped_files[i], chain_file=backwards_chain_file,
            mapped_file=backmapped_files[i], unmapped_file=unbackmapped_files[i]))

        # filter out intervals that does not map both ways
        gwf.target_from_template('{}_filter_{}'.format(slurm_tag, i),
            bed_difference(bed_file1=mapped_files[i], bed_file2=unbackmapped_files[i],
            output_file=filtered_files[i]))

    # filter out intervals that does not map both ways
    gwf.target_from_template('{}_merge_and_split'.format(slurm_tag),
        bed_merge_and_split(input_files=filtered_files, output_files=lifted_files))

    return lifted_files

    # # lift over each site individually 
    #sites_files = [liftover_dir / x.with_suffix('.sites').name for x in intervals_files]
    # for i, chrom in enumerate(sites_files):
    #   gwf.target('intr_to_sts_{}'.format(i)) << intervals_to_sites(intervals_file=intervals_files[i], sites_file=sites_files[i])
    #   gwf.target('lift_{}'.format(i)) << liftover(bed_file=sites_files[i], chain_file=forwards_chain_file,
    #       mapped_file=mapped_files[i], unmapped_file=unmapped_files[i])
    #   # lift back to orginal coordinates to ensure one to one correspondence
    #   gwf.target('liftback_{}'.format(i)) << liftover(bed_file=mapped_files[i], chain_file=backwards_chain_file,
    #       mapped_file=backmapped_files[i], unmapped_file=unbackmapped_files[i])
    #   # filter out intervals that does not map both ways
    #   gwf.target('liftbackfilter_{}'.format(i)) << bed_difference(bed_file1=mapped_files[i], bed_file2=unbackmapped_files[i],
    #       output_file=lifted_files[i])

    # # check that there are no overlapping (double) sites
    # uniq -cd <output_file>


def compute_interference_score(sites_files, annotation_files, map_files, scores_dir, map_scale, rate_scale):
    """
    Computes inteference scores for a specified set of
    sites given positions for interfering sites and a 
    recombination map
    """

    if not scores_dir.exists():
        os.makedirs(str(scores_dir))

    score_files = [scores_dir / x.with_suffix('.scores').name for x in sites_files]

    # input files grouped by chromosome
    file_tuples = zip(sorted(sites_files, key=lambda x: x.name),
                      sorted(annotation_files, key=lambda x: x.name),
                      sorted(map_files, key=lambda x: x.name),
                      sorted(score_files, key=lambda x: x.name))

    for i, (sites_file, annot_file, map_file, score_file) in enumerate(file_tuples):

        # all names are chr[chrom].bed and should abviously match up...
        assert sites_file.name == annot_file.name
        assert sites_file.name == map_file.name

        gwf.target_from_template('ic_interf_{}_{}'.format(species, i),
            interference_score(sites_file=str(sites_file),
             annot_file=str(annot_file), map_file=str(map_file), score_file=str(score_file),
             map_scale=map_scale, rate_scale=rate_scale))

    return score_files


#################################################################################
# Workflow
#################################################################################

chains_dir = Path('/project/hri/faststorage/data/chain_files')
human_chromosomes = ['chr{}'.format(x) for x in list(range(1,23)) + ['X']]
ape_chromosomes = ['chr1', 'chr2a', 'chr2b'] + ['chr{}'.format(x) for x in list(range(3,23)) + ['X']]
steps_root_dir = 'interference_score'

#########################################################
# Interference score for chimp
#########################################################

# species
species = 'chimp'
annotation_tag = 'gerp'

# input data
vcf_file = Path('/project/hri/faststorage/people/dce/results/VCF/Callable/chimp.callable.vcf')
annotation_files = [Path(x) for x in ['/project/hri/faststorage/data/SNPeff/db/gerp/GERP_hg19.bed']]
orig_map_files = list()
for p in Path('/project/hri/faststorage/data/great-ape-recombination-master/chimp_map').glob('*.hg18.txt.filtered'):
    if not p.name == 'chr2-map.hg18.txt.filtered':
        orig_map_files.append(p)


split_vcf_files = split_file(vcf_file, Path(os.getcwd(), 'steps', steps_root_dir, species, 'split_sites_files'))

# create bed files from vcf files
sites_files = format_vcf_to_bed(split_vcf_files, species,
        steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, 'converted_vcf_files'))

# create bed files from stevison's recombination map files
map_files = format_orig_map_files_to_bed(orig_map_files,
        steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, 'converted_map_files'))

# lift sites from hg19 to species coordinates
sites_files_lifted_to_species = reciprocal_liftover(sites_files, 
    forwards_chain_file=chains_dir/'hg19ToPanTro2.over.chain', 
    backwards_chain_file=chains_dir/'panTro2ToHg19.over.chain',
    slurm_tag='ic_siteslift_{}'.format(species),
    steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, 'sites_species_coord'),
    target_chromosomes=ape_chromosomes)

# lift annotation from hg19 to species coordinates
annotation_files_lifted_to_species = reciprocal_liftover(annotation_files,
    forwards_chain_file=chains_dir/'hg19ToPanTro2.over.chain', 
    backwards_chain_file=chains_dir/'panTro2ToHg19.over.chain',
    slurm_tag='ic_annotlift_{}'.format(species),
    steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, annotation_tag, 'annotation_species_coord'),
    target_chromosomes=ape_chromosomes)

# # lift recombination map from gorgor3 to hg19 coordinates
# map_files_lifted_to_hg19 = reciprocal_liftover(map_files,
#     forwards_chain_file=chains_dir/'gorGor3ToHg19.over.chain', 
#     backwards_chain_file=chains_dir/'hg19ToGorGor3.over.chain',
#     slurm_tag='ic_maplift1_{}'.format(species),
#     steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, 'map_hg19_coord'),
#     target_chromosomes=human_chromosomes)

# # lift recombination map from hg19 to species coordinates
# map_files_lifted_to_species = reciprocal_liftover(map_files_lifted_to_hg19,
#     forwards_chain_file=chains_dir/'hg19ToPanTro2.over.chain', 
#     backwards_chain_file=chains_dir/'panTro2ToHg19.over.chain',
#     slurm_tag='ic_maplift2_{}'.format(species),
#     steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, 'map_species_coord'),
#     target_chromosomes=ape_chromosomes)

# # sanity check (there should be one file for each chromosome)
# assert len(sites_files_lifted_to_species) == 24 # 1, 2a, 2b, 3, 4, ... 22, X
# assert len(sites_files_lifted_to_species) == len(annotation_files_lifted_to_species)
# assert len(sites_files_lifted_to_species) == len(map_files_lifted_to_species)


# compute interference scores
scores_files_species_coord = compute_interference_score(sites_files_lifted_to_species, 
    annotation_files_lifted_to_species, map_files,
    scores_dir = Path(os.getcwd(), 'steps', steps_root_dir, species, annotation_tag, 'scores_species_coord'),
    map_scale=1000,
    rate_scale=1) # use rho, not recombination rate

                      # Maybe instead normalize all maps to total length 1 
                      # in computation of genetic map.

# # lift scores back to hg19
# reciprocal_liftover(scores_files_species_coord, 
#         forwards_chain_file=chains_dir/'panTro2ToHg19.over.chain', 
#         backwards_chain_file=chains_dir/'hg19ToPanTro2.over.chain',
#         slurm_tag='ic_scoreslift_{}'.format(species),
#         steps_dir=Path(os.getcwd(), 'steps', steps_root_dir, species, annotation_tag, 'scores_human_coord'),
#         target_chromosomes=human_chromosomes)


#########################################################
# Interference score for bonobo
#########################################################

# species
species = 'bonobo'

# and so on....
 

#########################################################
# Interference score for gorilla
#########################################################

# species
species = 'gorilla'

# and so on....
 

#########################################################
# Interference score for gorilla
#########################################################

# species
species = 'orang'

# and so on....

# no recombination map for orangs yet :/

if __name__ == "__main__":
    print(len(gwf.targets))
