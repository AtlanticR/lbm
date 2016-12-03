log_gaussian_offset = function(offset=0) {
  structure(list(
    linkfun = function(mu) log(mu + offset), 
    linkinv = function(eta) exp(eta) - offset,
    mu.eta = function(eta) NA, 
    valideta = function(eta) TRUE, 
    name = paste0("logexp(", offset, ")") ),
    class = "link-glm" )
}
# or to make your own
  # conker_family_new = function(offset=0) {
  #   structure(list(
  #     linkfun = function(mu) mu + offset,
  #     linkinv = function(eta) mu - offset,
  #     mu.eta = function(eta) NA,
  #     valideta = function(eta) TRUE,
  #     name = paste0("logexp(", offset, ")") ),
  #     class = "link-glm" )
  # }
