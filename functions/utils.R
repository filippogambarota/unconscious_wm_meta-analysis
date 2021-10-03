#### UTILS FUNCTIONS

# rm_path_and_ext ---------------------------------------------------------

# remove path and extensions from a string(s)

rm_path_and_ext <- function(strings){
  strings %>% 
    basename() %>% 
    tools::file_path_sans_ext() 
}

# load_models ----------------------------------------------------------

# Wrapper to load all models (or a subset) from a specific folder and
# returning a named list. path = "mod" is the default path. It becomes
# irrelevant when the mod_names argument is not NULL. mod_names argument
# can be used to load a subset of models providing an array of models names

load_models <- function(mod_names = NULL, path = "mod"){
  
  if(is.null(mod_names)){
    mod_names <- list.files(path, pattern = ".rds", full.names = T)
  }
  
  mods <- map(mod_names, readr::read_rds)
  
  names(mods) <- rm_path_and_ext(mod_names)
  
  return(mods)
}


# add_to_list -------------------------------------------------------------

add_to_list <- function(list, element, name){
  list_names <- c(names(list), name)
  index <- length(list) + 1
  list[[index]] <- element
  names(list) <- list_names
  return(list)
}


# get_all_pakages ---------------------------------------------------------

# these fuctions greps all used packages among all scripts. The input is a
# vector of folders to explore and a vector of file extensions to detect

get_all_pakages <- function(folder_list, exts){
  folder_list <- here(folder_list)
  folder_list %>% 
    get_relevant_file(exts) %>% 
    get_all_files_content() %>% 
    detect_packages() %>% 
    get_pkg_name() 
}

get_relevant_file <- function(folder_list, exts){
  out <- unlist(map(folder_list, function(x) list.files(x, full.names = TRUE)))
  pattern <- paste0(exts, collapse = "|")
  index <- stringr::str_detect(out, pattern)
  out[index]
}

get_all_files_content <- function(files){
  unlist(map(files, function(x) readLines(x, warn = F)))
}

detect_packages <- function(files_content){
  index <- str_detect(files_content, "library\\(")
  out <- files_content[index]
  unique(out)
}

get_pkg_name <- function(pkgs){
  out <- str_extract(pkgs, pattern = "\\((.*?)\\)")
  out <- str_remove_all(out, "\\(|\\)")
  return(out)
}