
from itertools import islice
from math import exp

# import sys, traceback, os, multiprocessing


# class Error(Exception):
#     pass
    
# class ThreadException(Error):

#     def __init__(self, exception, traceback):
#         self.exception = exception
#         self.traceback = traceback

#     def __repr__(self):
#         return 'Error:\n' + str(self.exception) + "\n", self.traceback

# class Consumer(multiprocessing.Process):
    
#     def __init__(self, task_queue, result_queue):
#         multiprocessing.Process.__init__(self)
#         self.task_queue = task_queue
#         self.result_queue = result_queue

#     def run(self):
#         proc_name = self.name
#         while True:
#             next_task = self.task_queue.get()
#             if next_task is None:
#                 # Poison pill means we should exit
#                 break
#             try:
#                 result = next_task()
#             except Exception as exe:
#                 # This is a hack: we simply add the traceback to the result queue
#                 # and then check there if the result is a traceback
#                 result = 'Traceback:\n'+ "".join(traceback.format_tb(sys.exc_info()[2]))
#             self.result_queue.put(result)
#         return

# class Task(object):
    
#     def __init__(self, f, *args, **kwargs):
#         self.f = f
#         self.args = args
#         self.kwargs = kwargs
#     def __call__(self):
#         return self.f(*self.args, **self.kwargs)
    
# class Repeats(object):
    
#     def __init__(self, cores, num_jobs, f, *args, **kwargs):
#         self.tasks = multiprocessing.Queue()
#         self.results = multiprocessing.Queue()
#         self.num_jobs = num_jobs
#         self.cores = cores
#         self.f = f
#         self.args = args
#         self.kwargs = kwargs
        
#         #num_consumers = multiprocessing.cpu_count()
#         self.consumers = [ Consumer(self.tasks, self.results)
#                       for i in range(cores) ]

#         for w in self.consumers:
#             w.start()
    
#         self._queue_jobs()
            
#     def _queue_jobs(self):
    
#         # Enqueue jobs
#         for i in range(self.num_jobs):
#             self.tasks.put(Task(self.f, *self.args, **self.kwargs))

#         # Add a poison pill for each consumer
#         for i in range(self.cores):
#             self.tasks.put(None)
  
#     def __iter__(self):
#         return self

# #    def next(self): # Python 3: def __next__(self)
#     def __next__(self): # Python 3: def __next__(self)

#         if not self.num_jobs:
#             for w in self.consumers:
#                 w.join()
#             raise StopIteration
#         else:
#             while self.num_jobs:
#                 result = self.results.get()
#                 if type(result) is str and result.startswith('Traceback:'):
#                     del self.tasks
#                     del self.results
#                     print(result, file=sys.stderr)
#                     sys.stderr.flush()
#                     os.kill(os.getpid(),9)
#                     raise ThreadException
                    
#                 self.num_jobs -= 1
#                 return result

# class Broadcast(Repeats):

#     def __init__(self, cores, f, *args, **kwargs):

#         len_args, len_kwargs = 0, 0
#         if args:
#             len_args = len(args[0])
#             assert all(len(x) == len_args for x in args)
#         if kwargs:
#             vals = list(kwargs.values())
#             len_kwargs = len(vals[0])   # CHANGED THIS
#             assert all(len(x) == len_kwargs for x in vals)
#         if len_args and len_kwargs:
#             assert len_args == len_kwargs
#         assert len_args or len_kwargs, "Use the  class instead"

#         Repeats.__init__(self, cores, max(len_args, len_kwargs), f, *args, **kwargs)

#     # def __init__(self, cores, f, *args):

#     #     tups = list(zip(*args))

#     #     Repeats.__init__(self, cores, len(tups), f, *tups)

    
#     def _queue_jobs(self):
        
#         arg_list = list(zip(*self.args))

#         kwarg_list = list()
#         for i in range(len(list(self.kwargs.values())[0])):
#         # for i in range(len(self.kwargs)):
#             d = dict()
#             for k, v in self.kwargs.items():
#                 d[k] = v[i]
#             kwarg_list.append(d.copy())        
            
#         # Enqueue jobs
#         pairs = [[[], {}] for x in range(max(len(arg_list), len(kwarg_list)))]
#         for i in range(len(arg_list)):
#             pairs[i][0] = arg_list[i]
#         for i in range(len(kwarg_list)):
#             pairs[i][1] = kwarg_list[i]

#         assert len(pairs) == self.num_jobs, (len(pairs), self.num_jobs)
            
#         for args, kwargs in pairs:
#             self.tasks.put(Task(self.f, *args, **kwargs))
        
#         # Add a poison pill for each consumer
#         for i in range(self.cores):
#             self.tasks.put(None)
    
# ###########################################
    



def compute_scores(positions, stats, start_idx=None, end_idx=None):
    if start_idx is None: start_idx = 0
    if end_idx is None: end_idx = len(positions)
    for i, pos in enumerate(islice(positions, start_idx, end_idx)):
        scores = list()
        for f in stats:
            idx = i + start_idx
            downstream = sum(f(pos-p) for p in islice(positions, idx))
            upstream = sum(f(p - pos) for p in islice(positions, idx+1, len(positions)))
            scores.append(downstream + upstream)
        yield pos, scores

def inverse(d):
    return 1/d

def exponential(d):
    return exp(-d)

positions = list(range(100))

for pos, score in compute_scores(positions, stats=[exponential, inverse], start_idx=90, end_idx=95):
    print(pos, score)

# #for s in Broadcast(4, compute_scores, stats=[exponential, inverse], [4, 8, 4, 8], [14, 14, 14, 14]):
# for s in Broadcast(4, compute_scores, stats=[[exponential, inverse]], start_idx=[90], end_idx=[95]):
#     print(s)


    
# def square(number, constant):
#     return number**2 + constant

# for s in Broadcast(4, square, [4, 8, 4, 8], [14, 14, 14, 14]):
#     print(s)
# print()
# for s in Broadcast(1, square, number=[4, 8], constant=[14, 14]):
#     print(s)
# print()
# for s in Repeats(1, 4, square, 4, 14):
#     print(s)
# print()
# for s in Repeats(1, 2, square, number=4, constant=14):
#     print(s)
# print()


