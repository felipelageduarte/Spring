sourceCpp("cpp_source/gridSearch.cpp")


gridSearch <- function(series_obj, f, params){
    tic('gridSearch')
    gsResult = gridSearch_cpp(series_obj, f, params);
    #gsResult = gridSearch_r(series_obj, f, params);
    toc()
    return(gsResult)
}

gridSearch_r <- function(series_obj, f, params){
    series   = series_obj$series
    obsDet   = series_obj$det.series
    obsSto   = series_obj$sto.series
    gsResult = list()
    gsResult[["series_obj"]] = series_obj;
    gs = list()
    
    for (i in 1:nrow(params)){
        tic()
        det = f(series, params[i,]);
        t = toc(quiet = T)
        l = list();
        l[["series"]] = series;
        l[["params"]] = params[i,];
        l[["f"]]      = f;
        l[["det"]]    = det;
        l[["sto"]]    = series - det;
        l[["rt"]]     = t$toc - t$tic;
        gs[[i]] = l;
    }
    
    gsResult[["gs"]] = gs;
    gsResult
}

gridSearch_null_params <- function(series_obj, f){
    series   = series_obj$series
    obsDet   = series_obj$det.series
    obsSto   = series_obj$sto.series
    gsResult = list()
    gsResult[["series_obj"]] = series_obj;
    gs = list()
    
    tic()
    det = f(series, NULL);
    t = toc(quiet = T)
    l = list();
    l[["series"]] = series;
    l[["params"]] = NULL;
    l[["f"]]      = f;
    l[["det"]]    = det;
    l[["sto"]]    = series - det;
    l[["rt"]]     = t$toc - t$tic;
    gs[[i]] = l;

    gsResult[["gs"]] = gs;
    gsResult
}

runPipeline <- function(seriesList, f, params, folder, tech_name, run_r = FALSE){
    
    if(run_r == FALSE){
        sourceCpp("cpp_source/gridSearch.cpp")
    } else {
        cat("> Running in R mode.\n")
    }
    
    tic('pipeline')
    pb = txtProgressBar(1, length(seriesList), style=3, title=tech_name)
    for(i in 1:length(seriesList)){
        gs = list()
        if(is.null(params)){
            gs = gridSearch_null_params(seriesList[[i]], f)
        } else {
            if(run_r == TRUE){
                gs = gridSearch_r(seriesList[[i]], f, params)
            } else {
                gs = gridSearch_cpp(seriesList[[i]], f, params)
            }
        }
        save(gs, file=paste(folder,'/',tech_name, '_', sprintf('%04d', i), '.RData', sep=''))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    toc()
}