hypergeometric <- function(n1, m1, N, psi){    
  ##
  ## R program for computing the mean, variance, density, cumulative
  ## distribution and generating random deviates.
  ##
  ## Based on
  ## Liao and Rosen (2001): Fast and Stable Algorithms for Computing
  ## and Sampling from the Noncentral Hypergeometric Distribution,
  ## The American Statistician, 55, 236-369. 
  ##
  ## this is how to use the function
  ##
  ## n1 <- 100
  ## n2 <- 100
  ## m1 <- 100
  ## N <- n1+n2
  ## odds.ratio <- 3;
  ## obj <- hypergeometric(n1, m1, N, odds.ratio)
  ## obj$mean()
  ## obj$var()
  ## obj$d(40)
  ## obj$p(40)
  ## obj$r()
  
  n2 <- N - n1;
  
  if(n1<0 | n2<0 | m1<0 | m1>N | psi<=0)
    stop("wrong argument in hypergeometric");
  
  mode.compute <- function(){
    a <- psi - 1;
    b <- -( (m1+n1+2)*psi + n2-m1 ) ;    
    c <- psi*(n1+1)*(m1+1);
    q <- b + sign(b)*sqrt(b*b-4*a*c);
    q <- -q/2;
    
    mode <- trunc(c/q); 
    if(uu>=mode && mode>=ll) return(mode)
    else return( trunc(q/a) );      
  }       
  
  r.function <- function(i) (n1-i+1)*(m1-i+1)/i/(n2-m1+i)*psi;
  
  ##
  mean <- function() sum( prob[(ll:uu)+shift]*(ll:uu) ); 
  
  var <-  function() sum( prob[(ll:uu)+shift]*(ll:uu)^2 ) - mean()^2;          
  
  d <- function(x) return(prob[x + shift]);
  
  p <- function(x, lower.tail=TRUE){   
    if(lower.tail) return( sum(prob[ll:(x+shift)]) )
    else return( sum( prob[(x+shift):uu] ) );
  }
  
  ##
  
  sample.low.to.high <- function(lower.end, ran){ 
    for(i in lower.end:uu){                                
      if(ran <= prob[i+shift]) return(i);
      ran <- ran - prob[i+shift];
    }                                
  }
  
  sample.high.to.low <- function(upper.end, ran){           
    for(i in upper.end:ll){                              
      if(ran <= prob[i+shift]) return(i);
      ran <- ran - prob[i+shift];
    } 
  }  
  
  
  r <- function(){
    ran <- runif(1); 
    
    if(mode==ll) return( sample.low.to.high(ll, ran) );            
    if(mode==uu) return( sample.high.to.low(uu, ran) );                                         
    
    if(ran < prob[mode+shift]) return(mode);             
    ran <- ran - prob[mode+shift];
    
    lower <- mode - 1;                                                                            
    upper <- mode + 1;
    
    repeat{                                     
      if(prob[upper + shift] >= prob[lower + shift]){              
        if(ran < prob[upper+shift]) return(upper);
        ran <- ran - prob[upper+shift];
        if(upper==uu) return( sample.high.to.low(lower, ran) );
        upper <- upper + 1;                            
      }
      
      else{
        if(ran < prob[lower+shift]) return(lower);
        ran <- ran - prob[lower+shift];
        if(lower==ll) return( sample.low.to.high(upper, ran) );
        lower <- lower - 1;                   
      }      
    } 
  }
  
  ##
  
  ll <- max(0, m1-n2);
  uu <- min(n1, m1);              
  mode <- mode.compute();  
  
  prob <- array( 1, uu-ll+1 );
  
  shift <- 1-ll; 
  if(mode<uu) #note the shift of location
    {  
      r1 <- r.function( (mode+1):uu );       
      prob[ (mode+1 + shift):(uu + shift) ] <- cumprod(r1);       
    }
  
  if(mode>ll){
    r1 <- 1/r.function( mode:(ll+1) );
    prob[ (mode-1 + shift):(ll + shift) ] <- cumprod(r1);
  }
  
  prob <- prob/sum(prob); 
  
  return(list(mean=mean, var=var, d=d, p=p, r=r));
}
