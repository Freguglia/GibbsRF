complete_cMat = function(cMat){
  as.matrix(rbind(cMat,-cMat,deparse.level = 0)) %>% return
}

complete_vMat = function(vMat){
  A = ncol(vMat)
  (rbind(vMat,vMat[,A:1],deparse.level=0)) %>% return
}

gModel2vec = function(gModel,type){
  G = gModel$G
  C = gModel$cMat %>% nrow
  vMat = gModel$vMat
  V = gModel$V
  if(type=="general"){
    npars = G + C*(2*G)
    params = numeric(npars)
    if(!is.null(V) && !is.null(vMat)){
      params[1:G] = V[2:(G+1)]
      end = G
      for(i in 1:C){
        params[(end+1):(end+2*G)] = vMat[i,-(G+1)]
        end = end + 2*G
      }
    }
    return(params)
  }

  if(type=="symetric"){
    npars = G + C*(G)
    params = numeric(npars)
    if(!is.null(V) && !is.null(vMat)){
      end = G
      params[1:G] = V[2:(G+1)]
      for(i in 1:C){
        params[(end+1):(end+G)] = vMat[i,(G+2):(2*G + 1)]
        end = end+G
      }
    }
    return(params)
  }

  if(type=="equal"){
    npars = G + C
    params = numeric(npars)
    if(!is.null(V) && !is.null(vMat)){
      params[1:G] = V[2:(G+1)]
      params[(G+1):(G+C)] = vMat[,(G+2)]
    }
    return(params)
  }
}

vec2gModel = function(vec,cMat,G,type){
  V = numeric(G+1)
  V[2:(G+1)] = vec[1:G]
  V[1] = -sum(V)

  C = nrow(cMat)
  vMat = matrix(0,nrow=C,ncol=2*G+1)
  end = G
  for(i in 1:C){
    if(type=="general"){
      vMat[i,(1:G)] = vec[(end+1):(end+G)]
      vMat[i,((G+2):(2*G+1))] = vec[(end+G+1):(end+2*G)]
      vMat[i,(G+1)] = -sum(vMat[i,])
      end = end+2*G
    }
    if(type=="symetric"){
      vMat[i,(G+2):(2*G+1)] = vec[(end+1):(end+G)]
      vMat[i,1:G] = vec[(end+G):(end+1)]
      vMat[i,(G+1)] = -sum(vMat[i,(G+1):(2*G+1)])
      end = end+G
    }
    if(type=="equal"){
      vMat[i,-(G+1)] = vec[(end+1)]
      vMat[i,(G+1)] = -vec[(end+1)]
      end = end+1
    }
  }


  return(GibbsModel(G,cMat,V,vMat))
}
