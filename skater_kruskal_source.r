#HOW TO CALCULATE THE DISTANCE BETWEEN TWO NODES
mydist<-function(x1,x2){
	return((x1-x2)^2)
}
#RANDOM SPANING TREE BASED ON WEIGHTS (KRUSKAL)
rndtree<-function (nbw,wgt_1) 
{
	n <- length(nbw[[2]])
	
	from=c()
	to=c()
	wgt = c()
	
	finalfrom=c()
	finalto=c()
	finalwgt = c()
	
	conected = rep(0,n)
	
	#CONSTRUCT THE GRAPH CONSIDERING THE WEIGHTS(wgt_1)
	for(i in 1:n){
		nbw_i= nbw$neighbours[[i]]
		nbw_i_wgt= nbw$weights[[i]]
		
		from = c(from,rep(i,sum(nbw_i>i)))
		to = c(to,nbw_i[nbw_i>i])
		
	}
	
	for(i in 1:length(from)){
		wgt = c(wgt,mydist(wgt_1[from[i]],wgt_1[to[i]]))
	}
	
	
	#######################################################
	
	
	#ADD A EDGE IN EACH STEP
	while(sum(conected==1)!=length(conected) & (length(from)>0) ){	
		
		
		ac_wgt = 1/wgt #Higher weight, Lower probability
		ac_wgt = ac_wgt / sum(ac_wgt)
		ac_wgt = cumsum(ac_wgt)
		
		x = runif(1,0, max(ac_wgt))
		idx = (1:length(ac_wgt))[ac_wgt==(ac_wgt[ac_wgt>x])[1]]
		idx = idx[1] #Choosing randomily the edge to put in the Random Tree
	
		from_idx = from[idx]
		to_idx = to[idx]
		wgt_idx = wgt[idx]
	
		#Is this one of those nodes connected to the graph?
		if(conected[from_idx]==0 & conected[to_idx]==0){
			#Create a new tree 
			newgroup = max(conected)+1
			conected[from_idx] = newgroup
			conected[to_idx] = newgroup
		
			finalfrom = c(finalfrom,from_idx)
			finalto = c(finalto,to_idx)
			finalwgt = c(finalwgt, wgt_idx)
		
		}else{
			if(conected[from_idx]!=conected[to_idx]){
				#Connect trees
				maxft = max(conected[from_idx],conected[to_idx])
				minft = min(conected[from_idx],conected[to_idx])
				if(minft == 0) minft = maxft
			
				conected[from_idx] = minft
				conected[to_idx] = minft
			
				conected[conected == maxft] = minft
			
				finalfrom = c(finalfrom,from_idx)
				finalto = c(finalto,to_idx)
				finalwgt = c(finalwgt, wgt_idx)
			}
	
		}	
	
		#Take out the edge
		from = from[-idx]
		to = to[-idx]
		wgt = wgt[-idx]
		
		
	
	}
	
	rntree = matrix(c(finalfrom,finalto,finalwgt),ncol = 3)
	attr(rntree, "class") <- c("mst", "matrix")
    rntree
	
}

#KRUSKAL ALGORITHM
#Algorithm equivalent to rndtree (Changing just the way to select a edge)
kruskal<-function (nbw,wgt_1) 
{
	n <- length(nbw[[2]])
	
	from=c()
	to=c()
	wgt = c()
	
	finalfrom=c()
	finalto=c()
	finalwgt = c()
	
	conected = rep(0,n)
	
	for(i in 1:n){
		nbw_i= nbw$neighbours[[i]]
		nbw_i_wgt= nbw$weights[[i]]
		
		from = c(from,rep(i,sum(nbw_i>i)))
		to = c(to,nbw_i[nbw_i>i])
		
	}
	
	for(i in 1:length(from)){
		wgt = c(wgt,mydist(wgt_1[from[i]],wgt_1[to[i]]))
	}
	
	
	
	while(sum(conected==1)!=length(conected) & (length(from)>0) ){	
		

		
		ac_wgt  = wgt 
		idx = (1:length(ac_wgt))[(ac_wgt==min(ac_wgt))]
		idx = idx[1] #Selecting the edge with the minimum weight
	
		from_idx = from[idx]
		to_idx = to[idx]
		wgt_idx = wgt[idx]
	
		if(conected[from_idx]==0 & conected[to_idx]==0){
			newgroup = max(conected)+1
			conected[from_idx] = newgroup
			conected[to_idx] = newgroup
		
			finalfrom = c(finalfrom,from_idx)
			finalto = c(finalto,to_idx)
			finalwgt = c(finalwgt, wgt_idx)
		
		}else{
			if(conected[from_idx]!=conected[to_idx]){
				maxft = max(conected[from_idx],conected[to_idx])
				minft = min(conected[from_idx],conected[to_idx])
				if(minft == 0) minft = maxft
			
				conected[from_idx] = minft
				conected[to_idx] = minft
			
				conected[conected == maxft] = minft
			
				finalfrom = c(finalfrom,from_idx)
				finalto = c(finalto,to_idx)
				finalwgt = c(finalwgt, wgt_idx)
			}
	
		}	
	
		
		from = from[-idx]
		to = to[-idx]
		wgt = wgt[-idx]
		
		
	
	}
	
	rntree = matrix(c(finalfrom,finalto,finalwgt),ncol = 3)
	attr(rntree, "class") <- c("mst", "matrix")
    rntree
	
}

#FUNCTION TO REMOVE A EDGE FOR A GRAPH AND SEPARATE THE TREES AS CLUSTERS
newclusters<-function(clusters,edge,from,to){
	
	
	newcluster = max(clusters) + 1
	f1 = edge[1]
	t1 = edge[2]
	
	pcluster = clusters[f1]

	
	tos = c(f1,t1)
	vis = c(0,1)
	
	
	
	
	while(sum(vis)>0){
		at = ((1:length(vis))[(vis==1)])[1]
		vis[at] = 0
		
		node = tos[at]
		
		if(clusters[node]==pcluster){
			clusters[node] = newcluster
			ndx = to[(1:length(from))[(from==node)]]
			ndx = c(ndx,from[(1:length(to))[(to==node)]])
			ndx = unique(ndx)
			vis=c(vis,rep(1,length(setdiff(ndx,tos))))
			tos = c(tos,setdiff(ndx,tos))
		}
		
		
	}
	
	
	return(clusters)
	
}

#SKATER - PRUNE ALGORITHM
mySkater<-function(tree,stat,prune){
	clusters = rep(1,length(stat))
	
	from = tree[,1]
	to = tree[,2]
	ssws = c(var(stat))
	
	#Each iteration prunes the tree by removing an edge
	for(k in 1:prune){
		vars = c()
		#Try each possible an edge
		for(i in 1:length(from)){

			nclusters =  newclusters(clusters,c(from[i],to[i]),from,to)
			n = max(nclusters)
			vartot = 0;
			for(j in 1:n){
				if(length(stat[nclusters==j])>1)
					vartot = vartot + var(stat[nclusters==j])
			}
			vars = c(vars,vartot)
		}
		
		idx = (1:length(from))[(vars == min(vars))]
		idx = idx[1] #Selects the edge that holds the lowest variance
		
		ssws=c(ssws,vars[idx])
		
		clusters = newclusters(clusters,c(from[idx],to[idx]),from,to)
		
		from=from[-idx]
		to=to[-idx]
		
	}
	
	res = list(
		regions = clusters,
		ssw = ssws
	)
	
	
	
	return(res)
	
}

#RANDOM SKATER (PRUNE) BASED ON WEIGHTS
#Algorithm equivalent to myskater (Changing just the way to select a edge)
mySkater_rnd<-function(tree,stat,prune){
	clusters = rep(1,length(stat))
	
	from = tree[,1]
	to = tree[,2]
	ssws = c(var(stat))
	
	var2 = c()
	
	for(k in 1:prune){
		vars = c()
		for(i in 1:length(from)){

			nclusters =  newclusters(clusters,c(from[i],to[i]),from,to)
			n = max(nclusters)
			vartot = 0;
			for(j in 1:n){
				if(length(stat[nclusters==j])>1)
					vartot = vartot + var(stat[nclusters==j])
			}
			vars = c(vars,vartot)
		}		
		
		ac_wgt = 1/vars  #Higher variance, Lower probability
		ac_wgt = ac_wgt / sum(ac_wgt)
		ac_wgt = cumsum(ac_wgt)
		
		x = runif(1,0, max(ac_wgt))
		idx = (1:length(ac_wgt))[ac_wgt==(ac_wgt[ac_wgt>x])[1]]
		idx = idx[1] #Choosing randomily the edge to remove of the Spanning Tree
		
	
		
		ssws=c(ssws,vars[idx])
		
		clusters = newclusters(clusters,c(from[idx],to[idx]),from,to)
		
		from=from[-idx]
		to=to[-idx]
		
	}
	
	res = list(
		regions = clusters,
		ssw = ssws
	)
	
	
	
	return(res)
	
}

