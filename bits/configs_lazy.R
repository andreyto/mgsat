test1 <- function() {
  library(pryr)
  a = within(list(), {
    b %<d-% 1
    c = b
    d %<d-% within(list(), {
      e %<d-% (b*2)
    })
  })
  
  p = within(a,{
    b %<d-% 2
  })
  print(p$d)
}

test2 <- function() {
  library(proto)
  a = proto(b=1)
  a$d = a$b*2
  p = a$proto(b=2)
  print(p$d)
}

test3 <- function() {
  library(proto)
  a = new.env()
  with(a, {
    b = 1
    c = b
    d = new.env()
    makeActiveBinding("e", function() (b*2),d)
  })
  p = a
  #p = a
  with(p,{
    b = 2
  })
  print(p$d$e)
}

test4 <- function() {
  library(proto)
  a = proto()
  with(a, {
    b = 1
    c = b
    d = proto()
    makeActiveBinding("e", function() (b*2),d)
  })
  p = a$proto()
  with(p,{
    b = 2
  })
  print(p$d$e)
}

#test5 <- function() {
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

AB <- function(env = parent.frame(),...) {
  this_call = match.call()
  #env = parent.frame()
  for(i in 2:length(this_call)) {
    var_name = names(this_call)[[i]]
    #cat("var_name=",var_name)
    if (exists(var_name, envir = env, inherits = FALSE)) {
      rm(list = var_name, envir = env)
    }
    makeActiveBinding(var_name, make_function(alist(x = ), this_call[[i]], env),env)
  }
}

EV <- function(.with=NULL,
               .update=T,
               .under=parent.frame(),...) {
  #cat(str(names(sys.call())))
  this_call = match.call()
  if(is.null(.with)) {
    env = new.env(parent=.under)
  }
  else {
    env = .with
  }
  for(i in 2:length(this_call)) {
    var_name = names(this_call)[[i]]
    if(substring(var_name,1,1)!=".") {
      var_name_exists = exists(var_name, envir = env, inherits = FALSE)
      var_rhs = this_call[[i]]
      is_rhs_ev = (var_rhs[[1]]=="EV")
      if(is_rhs_ev) {
        ## Missing .update arg means T (the default, but default is not set in
        ## the unevaluated expression yet)
        update_ev = (is.null(var_rhs$.update) || eval(var_rhs$.update,env))
        if(!var_name_exists) {
          update_ev = F
        }
        if(update_ev) {
          var_rhs[[".with"]] = eval(parse(text=var_name),env)
        }
      }
      func_rhs = make_function(alist(x = ), var_rhs, env)
      #if (exists(var_name, envir = env, inherits = FALSE)) {
      #  rm(list = var_name, envir = env)
      #}
      if(is_rhs_ev) {
        env[[var_name]] = func_rhs()
      }
      else {
        makeActiveBinding(var_name, func_rhs, env)
      }
    }
  }
  env
}

copy_env <- function(x,deep=F,
                     parent=parent.env(x),
                     deep.update.parent=T,
                     envir=NULL,
                     ...) {
  y = list2env(as.list(x, all.names=TRUE),
               parent = parent,envir=envir,...)
  if(deep) {
  for (name in names(y)) {
    if(is.environment(y[[name]])) {
      y_sub = y[[name]]
      y_sub_par = parent.env(y_sub)
      if(deep.update.parent && identical(y_sub_par,x)) {
        y_sub_par = y
      }
      y[[name]] = copy_env(y_sub,deep = deep,parent = y_sub_par,
                           deep.update.parent = deep.update.parent,
                           ...)
      ##TODO: update frame of functions and active bindings with environment()<-
    }
  }
  }
  return (y)
}

test6 <- function() {
  evlist <- function(expr,ev.par=parent.frame()) {
    expr = substitute(expr)
    e = new.env(parent=ev.par)
    with(e,{
      `<-` <- function(x,y) makeActiveBinding(substitute(x), make_function(alist(x = ), substitute(y), parent.frame()),parent.frame())
    })
    eval(expr,e)
    e
  }
  a = evlist({
    b <- 1
  })
  s = evlist({
    v <- 1 + b
    y <- if(b==1) (v + 1) else (v + 3)
    c = b
    AB(w = b*2)
    hidden = evlist({
      AB(z = w*3)
      w = 1
    })
  },a)
  with(a, 
       {b <- 2}
  )
}

test7 <- function() {
EV(s=4,t=EV(x=s*2))$t$x
x = EV(s=1:4,
       f=function(x) sprintf("x is %s",x),
       t=EV(x=s*2,
            z=EV(k=t$x*4)
       )
)
print(x$t$z$k)
if(T) {
  y = EV(.with = x,
         s=1:16,
         f=function(x) sprintf("y is %s",x),
         t=EV(l=10))
  print(y$f(y$t$z$k))
}
}

test8 <- function() {
  f2 <- function() {
    v1 <- 1
    v2 <- 2
    #environment(f1) <<- environment()
    return (function() {})
  }
  print(as.list(environment(f2())))
  return (f2())
}

