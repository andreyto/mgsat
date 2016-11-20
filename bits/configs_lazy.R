#library(pryr)
to_env <- function (x, quiet = FALSE) 
{
  if (is.environment(x)) {
    x
  }
  else if (is.list(x)) {
    list2env(x)
  }
  else if (is.function(x)) {
    environment(x)
  }
  else if (length(x) == 1 && is.character(x)) {
    if (!quiet) 
      message("Using environment ", x)
    as.environment(x)
  }
  else if (length(x) == 1 && is.numeric(x) && x > 0) {
    if (!quiet) 
      message("Using environment ", search()[x])
    as.environment(x)
  }
  else {
    stop("Input can not be coerced to an environment", call. = FALSE)
  }
}

all_named <- function (x) 
{
  if (length(x) == 0) 
    return(TRUE)
  !is.null(names(x)) && all(names(x) != "")
}

make_function <- function(args, body, env = parent.frame()) {
  args <- as.pairlist(args)
  #stopifnot(
  #  all_named(args),
  #  is.language(body))
  env <- to_env(env)
  
  eval(call("function", args, body), env)
}


rm_if_exists <- function(name,env) {
  for(nm in name) {
    if (exists(nm, envir = env, inherits = FALSE)) {
      rm(list = nm, envir = env)
    }
  }
} 

acon <- function(.with=NULL,
                 .update=T,
                 .under=parent.frame(),...) {
  ## How this is designed:
  ## The general idea is to overload the `=` assignment operator.
  ## Since it is impossible in R, we instead emulate it by 
  ## designing a special processing of the `=` when it is used
  ## to specify values for named arguments in a function call.
  ## this acon() constructor performs NSE of its argument list and
  ## returns a new environment object (S3 subclassed as `acon`).
  ## Any value expression that is itself a call to acon() is generated with a recursive call
  ## to acon(), but with the .under argument updated to point to the currently built
  ## environment (current acon object).
  ## Any other value expression is first evaluated anc checked if it is of type acon. If yes,
  ## the evaluated value is copied and placed under the current acon. If no, the name of
  ## the argument is made an active binding bound to a function that is constructed from the
  ## value expression.
  ## Thus, the returned result is a nested set of environments (acon objects), wit leaf values
  ## all being active bindings. When any of the active binding attributes is accessed, the
  ## corresponding expression will be evaluated, with names in the expression automatically being
  ## looked up the chain of nested environments until found.
  this_call = match.call()
  if(is.null(.with)) {
    env = new.env(parent=.under)
  }
  else {
    ## always make a deep copy of the .with environment, thus providing a "copy
    ## constructor" semantics
    env = copy_env(.with,deep = T,parent = .under)
  }
  class(env) <- append(class(env),"acon",0)
  #   print("DEBUG START")
  #   print(paste0(".parent= ",as.list(.under)," address_parent=",data.table::address(.under),
  #                " this_call=",paste(this_call,collapse = ",")," address_env=",data.table::address(env),
  #                " sys.calls=",paste(sys.calls(),collapse = "->")))
  #   print("DEBUG END")
  
  for(i in 2:length(this_call)) {
    var_name = names(this_call)[[i]]
    if(substring(var_name,1,1)!=".") {
      var_name_exists = exists(var_name, envir = env, inherits = FALSE)
      var_rhs = this_call[[i]]
      if(is.language(var_rhs) && (!is.symbol(var_rhs)) && var_rhs[[1]]=="acon") {
        ## Missing .update arg means T (the default, but default is not set in
        ## the unevaluated expression yet)
        update_ev = (is.null(var_rhs$.update) || eval(var_rhs$.update,env))
        ## if rhs is acon ctor expression, and the lhs already exists in the
        ## current object (inherited from .with), then the lhs is assumed to be
        ## the environment and is used as the .with argument to call the rhs ev ctor.
        if(!var_name_exists) {
          update_ev = F
        }
        if(update_ev) {
          var_rhs[[".with"]] = eval(parse(text=var_name),env)
        }
        var_rhs[[".under"]] = env
        ## simply evaluate the acon ctor call and assign
        val = eval(var_rhs,envir = env)
        rm_if_exists(var_name,env)
        env[[var_name]] = val
      }
      else {
        ## we have to eval in order to check the class of the resulting value
        val = eval(var_rhs,env)
        ## if rhs is acon object, it is copied and attached to the current object
        if(inherits(val,"acon")) {
          val = copy_env(val,deep = T,parent = env)
          rm_if_exists(var_name,env)
          env[[var_name]] = val
        }
        ## if not acon object or acon ctor expression, assign as active binding
        else {
          func_rhs = make_function(alist(), var_rhs, env) #.x = 
          rm_if_exists(var_name,env)
          makeActiveBinding(var_name, func_rhs, env)
        }
      }
    }
  }
  env
}

as.list.acon <- function(x) {
  as_list_nested_env(x)
}

as_list_nested_env <- function(x) {
  done_env_hash = new.env()
  
  as_list_nested_env_inner <- function(x) {
    stopifnot(is.environment(x) || is.list(x))
    y = list()
    x_names = names(x) #for list elements w/o names, returns ""
    for(i_el in seq_along(x_names)) {
      if(is.environment(x)) val = x[[x_names[[i_el]]]]
      else val = x[[i_el]]
      if(is.environment(val)) {
        val_key = data.table::address(val)
        if( is.null(done_env_hash[[val_key]]) ) {
          done_env_hash[[val_key]] = T
          val = as_list_nested_env_inner(val)
        }
      }
      else if(is.list(val)) {
        val = as_list_nested_env_inner(val)
      }
      y[[i_el]] = val
    }
    names(y) = names(x)
    return (y)
  }
  return (as_list_nested_env_inner(x))
}

ls_str_nested_env <- function(x) {
  ls.str(as_list_nested_env(x))
}

print_nested_env <- function(x,max.level=100,...) {
  print(ls_str_nested_env(x),max.level=max.level,...)
}

copy_env <- function(x,deep=F,
                     parent=parent.env(x),
                     func.update.parent=T,
                     deep.update.parent=T,
                     into_envir=NULL,
                     ...) {
  ## as.list copies active binding as functions; this allows us 
  ## recreating them as active bindings in the new environment
  y = list2env(as.list.environment(x, all.names=TRUE),
               parent = parent,envir=into_envir,...)
  for (name in names(y)) {
    if(deep) {
      
      if(is.environment(y[[name]])) {
        y_sub = y[[name]]
        y_sub_par = parent.env(y_sub)
        if(deep.update.parent && identical(y_sub_par,x)) {
          y_sub_par = y
        }
        y[[name]] = copy_env(y_sub,deep = deep,parent = y_sub_par,
                             func.update.parent = func.update.parent,
                             deep.update.parent = deep.update.parent,
                             ...)
      }
    }
    ## both plain functions and active bindings become functions in y
    if(is.function(y[[name]])) {
      y_sub = y[[name]]
      y_sub_par = environment(y_sub)
      ## update function environment if the original one points to the enclosing object
      if(func.update.parent && identical(y_sub_par,x)) {
        y_sub_par = y
        environment(y_sub) <- y_sub_par
      }
      ## additionally, if the original value was active binding, create active binding for the
      ## copied value
      if(bindingIsActive(name,x)) {
        rm(list = name, envir = y)
        makeActiveBinding(name,y_sub,y)
      }
      else {
        y[[name]] = y_sub
      }
    }
  }
  return (y)
}

test_acon_constructor <- function() {
  z = acon(v=4,w=acon(x=v*2))
  x = acon(s=1:4,
           f=function(x) sprintf("x is %s",x),
           t=acon(x=s*2,
                  z=acon(k=x*4)
           ),
           p = z
  )
  print(x$t$z$k)
  y = acon(x,
           s=1:16,
           f=function(x) sprintf("y is %s",x),
           t=acon(l=10))
  print(y$f(y$t$z$k))
  return (list(y=y,x=x))
}

test_copy_env <- function() {
  x = list2env(list(a=2,b="env_x"))
  x$y = list2env(list(c=2,d="env_y"),parent = x)
  x$y$z = list2env(list(e=2,f="env_z"),parent = x$y)
  print("x=")
  print_nested_env(x)
  a = copy_env(x,deep = T)
  a$a = 3
  a$y$c = 3
  a$y$z$e = 3
  print("a=")
  print_nested_env(a)
}

