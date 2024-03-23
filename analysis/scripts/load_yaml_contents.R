load_yaml_contents <- function(yaml_file_path) {
    config <- yaml::read_yaml(yaml_file_path)
    
    for (name in names(config)) {
        value <- config[[name]]
        
        # Check if the value starts with "glue::glue"
        if (class(value) == "character" && startsWith(value, 'glue::glue')) {
            # Extract the expression inside the glue call
            expr <- sub('glue::glue\\(["\'](.*)["\']\\)', "\\1", value)
            
            # Evaluate the expression within the glue function
            # and assign it globally
            evaluated_value <- glue::glue(expr)
            assign(name, evaluated_value, envir = .GlobalEnv)
        } else {
            assign(name, value, envir = .GlobalEnv)
        }
    }
}