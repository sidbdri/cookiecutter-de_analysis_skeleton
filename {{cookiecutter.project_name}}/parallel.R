NUM_CORES <- {{cookiecutter.number_total_threads}}

set_num_cores <- function(cores) {
  set_global(cores, "NUM_CORES")
}

run_jobs_with_shared_memory <- function(job_strings, 
                                        job_function='run_topgo',
                                        job_source = NULL,
                                        ...) {
 run_jobs(job_strings, job_function, job_source, "MulticoreParam", ...)
}

run_jobs_with_separate_memory <- function(job_strings, 
                                          job_function='run_topgo',
                                          job_source = NULL,
                                          ...) {
  run_jobs(job_strings, job_function, job_source, "SnowParam", ...)
}

run_jobs <- function(job_strings, 
                     job_function='run_topgo',
                     job_source = NULL,
                     parallelParam=c('SnowParam','MulticoreParam')[1],
                     nc=NUM_CORES, # might be better to pass the NUM_CORES
                     ...) {

  message(length(job_strings),' jobs received for ', job_function)
  print(job_strings)
  
  log_dir <- file.path('results/logs/R/BioParallel', job_function)
  message('logs can be found at ', log_dir)
  dir.create(log_dir, recursive = T, showWarnings = F)
  
  # nc <- NUM_CORES
  if (length(job_strings) < nc) {
    message('the number of cores is set to ', nc, ', but there are fewer jobs than cores; reduce the number of cores...')
    nc = length(job_strings)
  }
  
  if (nc > 1) {
    message('running bioparallel in parallel with <', parallelParam, '>. number of cores: ', nc)
    
    para <- call(parallelParam, workers = nc, stop.on.error = TRUE, jobname = job_function, progressbar = T,
                 log = TRUE, threshold = "DEBUG", logdir = log_dir, tasks = length(job_strings)) %>% eval
    print(para)
    
    if (parallelParam == 'SnowParam') {
      message('init workers...')
      bpstart(para)
      
      # we start the workers this way so we can reuse them to hopefully decrease some overhead of loading packages
      # need to load package, as in SOCK mode, workers are independent
      # if we are reusing a worker, we don't need to repeat the package load
      prepare_worker <- function(...) {
        source('load_packages.R')
        source('utility_functions.R')
        
        if (!is.null(job_source)) {
          source(job_source)
        }
        
        'success'
      }
      
      message('preloading packages and sourcing required functions in workers...')
      loadpackage <- tryCatch(
        {
          bplapply(c(1:nc), prepare_worker, BPPARAM = para)
        },
        error = function(e){ bpstop(para); e}
      )
      check_bpresult(attr(loadpackage,'result'), job_name = 'load package/source common_function in worker node')
    }
    
    if (parallelParam == 'MulticoreParam') {
      message('no need to prepare workers for MulticoreParam...')
    }
    
  } else {
    message('running bioparallel in non-parallel mode with SerialParam')
    para <- SerialParam(stop.on.error = TRUE, log = TRUE, threshold = "DEBUG", progressbar = T, logdir = log_dir)
    print(para)
  }
  
  message('running jobs in workers....')
  tictoc::tic(paste0(job_function,' run time:'))
  ret <-  tryCatch(
    {
      bplapply(job_strings,
               job_function,
               ...,
               BPPARAM = para
      ) %>% set_names(job_strings)
    },
    error = function(e) e, finally = bpstop(para)) # we need to stop the workers when error
  tictoc::toc()
  
  check_bpresult(attr(ret, 'result'), job_name = job_function)
  ret
}

check_bpresult <- function(res,job_name) {
  if (!all(bpok(res))) {
    # at least one worker node has error.
    message(job_name, 'parallel run failed')
    message('Error in worker node ', which(!bpok(res)) %>% paste0(., collapse = ','))
    message('Error log can be found at: ', paste0('results/logs/R/BioParallel/',job_name))
    sapply(which(!bpok(res)), function(i) {
      message('Worker ',i,' error message:')
      print(res[[i]])
      message('Worker ',i,' error trackback:')
      print(attr(res[[i]], 'traceback'))
    }) %>% invisible() 
    stop(job_name, 'parallel run failed')
  }
  message(job_name,' parallel run success')
}