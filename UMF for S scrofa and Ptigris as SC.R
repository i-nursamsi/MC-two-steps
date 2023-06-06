######### Ultra-Loop to Create UMFs ##########

## The unmarked frame object (UMF) lies at the heart of all hierarchical models in 'unmarked'

## A UMF has space for -
### 1) a detection history matrix, called "dh.mat" in the code below, from caps file

### 2) sampling unit-level covariates, which vary per sampling-unit (i.e. row)
## These are called "m" in the code below, from meta file

### 3) sampling occasion-level covariates, which vary per sampling-occasion (i.e. column)
## These are called 'obs3' in the code below, from caps file 


# Select which species you want
species = c("Sus_scrofa", "Panthera_tigris") # Bread and butter of the ECL. Also one common, one rare

# specify if you want detection (for occu models) matrix or count (for abundance model) matrix
type = "occu"
# type = "abundance"


### We will compress the detection history matrix in order to decrease the number of 
### non-detections (i.e. 0's) and increase detection probability. 
### The larger the sampling occasion window, the less certainty we have for 
### calculating detection probablity, which generates a larger SE. 


## Select sampling occasion window length, you can change it!
w = 5
## I am using a 5 day window b/c that's what I used in the co-abundance MS

## Select maximum number of days the matrix can be active (i.e. total duration)
## To satisfy model assumptions about population closure during sampling. 

dur = 200
## I am using a 200 day duration b/c that's what I used in the co-abundance MS


###### The Ultra-Loop

umf.list = list()

# t=1; j=1; s=1 #can use these variable numbers to test individual lines in the loop

for(t in 1:length(species)){ #run this loop for each species
  
  sp = (species)[t] # select a single species
  
  ### Specify which survey the species has been detected (at least once!)  
  ## Dont use landscape b/c it could introduce too many zeros for unmarked 
  survs = unique(caps$survey_id[caps$Species == sp])
  
  ## Select relevant surveys
  c = caps[caps$survey_id %in% survs,] #subset captures
  m = meta[meta$survey_id %in% survs,] #subset metadata
  
  
  ## standardize site covaraites to ensure variables are comparable across models later
  m.num<- m[,sapply(m, is.numeric)] 
  m.std<- decostand(m.num, method = "standardize", na.rm = TRUE)
  m.std = m.std[,colSums(is.na(m.std)) < nrow(m.std)] #remove any columns with no data
  m.char<- m[,sapply(m, is.character)]
  m<- data.frame(m.char, m.std)
  
  
  #
  ##
  ###
  #### Begin Creating Detection/Count History Matrix 
  ###
  ##
  #
  
  
  ## Outline the structure (i.e. no data) of the matrix and add col and row names
  mat = matrix(NA, 
               nrow = length(unique(c$cell_id_1km)), # number of rows are our sampling locations
               ncol = length(seq(from =1, to= max(c$seq))), # number of columns are our sampling occasions
               dimnames = list(as.character(unique(c$cell_id_1km)), # row names, then column names
                               seq(from =1, to= max(c$seq))))
  
  ## Determine when each sampling unit was active-
  for(j in 1:length(unique(c$cell_id_1km))){
    
    a= c[c$cell_id_1km == unique(c$cell_id_1km)[j],] #subset for a specific sampling unit out of all possible options
    
    indx = seq(from = 1, to = max(a$seq)) #determine the sequence it was operational for 
    
    mat[j,indx]=0 # at row J, across all sampling occasions, put a zero there
  }
  
  ## Fill in the matrix based on sampling unit and sampling occasion
  for(j in 1:length(unique(c$cell_id_1km))){ #repeat for each sampling unit
    
    su = unique(c$cell_id_1km)[j] #specify the sampling unit
    
    a = c[c$cell_id_1km == su & c$Species == sp,] #subset captures data for specific sampling unit and for specific species
    
    
    
    # Fill in matrix w/ count data 
    if(nrow(a)>0 & type == "abundance"){ #Bypass cameras without a detection and leave them as zero
      
      for(s in 1:length(a$Date)){ #repeat for each Date detected at the sampling unit
        
        d = a$Date[s] #specify the sampling date
        
        indx = a$seq[a$cell_id_1km == su & a$Date == d] #specify the matching date index
        
        mat[su,indx]= a$total_indiv_records[a$Date == d]
        #in the matrix where the row = sampling unit, and column = occasion, 
        #use the total counts of the specific sampling occasion
        
      } # end count fill
    } # end if-abundance statement
    
    
    
    # Fill in matrix w/ presence data
    if(nrow(a)>0 & type == "occu"){ #Bypass cameras without a detection and leave them as zero
      
      for(s in 1:length(a$Date)){ #repeat for each Date detected at the sampling unit
        
        d = a$Date[s] #specify the sampling date
        
        indx = a$seq[a$cell_id_1km == su & a$Date == d] #specify the matching date index
        
        mat[su,indx]= 1
        #in the matrix where the row = sampling unit, and column = occasion, 
        #mark the presence of the species with a 1
        
      } # end presence fill
    } # end if-occu statement
    
  }# end matrix filling loop 
  
  
  #
  ##
  ### Compress Detection/Count History Matrix into multi-day Sampling Occasions
  ##
  #
  
  
  ## If all sampling units are active for less than 200 days, use the maximum sequence length for a instead
  if(max(c$seq) < dur){ dur = max(c$seq)}
  
  
  ## Create a new and empty compressed matrix to fit sampling occasions
  dh.mat = matrix(nrow = nrow(mat), ncol = round(dur/w))
  
  
  ### Sampling occasion loop-
  for(u in 1:nrow(mat)){ #Repeat for each row in the matrix
    
    for(p in 1:round(dur/w)){ # Repeat for each sampling occasion window 
      
      # Outline the start dates of sampling occasion, using values provided in for-loop
      starts<-seq(from=1, to=dur, by=w)
      
      # Select a single start date
      l<-starts[p]
      
      
      if(type == "abundance"){
        
        # Make if-else statement 
        ifelse(all(mat[u,c(l:(l+w))]== "NA", na.rm=FALSE) == "TRUE", 
               #if all values in matrix @ row u across the sampling occasion window are all NA,
               dh.mat[u,p]<-NA, # then leave the sampling occasion as NA,
               dh.mat[u,p]<- sum(as.numeric(mat[u,c(l:(l+w))]), na.rm = TRUE)) 
        # But if FALSE, take the the sum of the detections
        
      } # End conditional for abundance
      
      if(type == "occu"){
        
        # Make if-else statement 
        ifelse(all(mat[u,c(l:(l+w))]== "NA", na.rm=FALSE) == "TRUE", 
               #if all values in matrix @ row u across the sampling occasion window are all NA,
               dh.mat[u,p] <- NA, # then leave the sampling occasion as NA,
               dh.mat[u,p] <- max(as.numeric(mat[u,c(l:(l+w))]), na.rm = TRUE))
        # But if FALSE, take the the maximum value (either zero or one)
        
      } # End conditional for occupancy
      
    } # End loop per sampling occasion
    
  } # End loop per row in matrix 
  
  
  #
  ##
  ### Format Observation covariates to match compressed matrix 
  ##
  #
  
  
  ## Generate observation covariates to match sampling occasions 
  obs = distinct(dplyr::select(c, cell_id_1km, seq, 
                               num_cams_active_at_date)) ## Come here and change to include more obs.covs if we get them!
  
  ## create empty obs dataframe
  obs.mat = (matrix(NA, 
                    nrow = length(unique(obs$cell_id_1km)), 
                    ncol = length(seq(from =1, to= max(obs$seq))),
                    dimnames = list(as.character(unique(obs$cell_id_1km)), # row names, then column names
                                    seq(from = 1, to= max(obs$seq)))))
  
  for(u in 1:length(unique(obs$cell_id_1km))){ #repeat for each sampling unit
    
    # Select a single sampling unit (i.e. row)
    su = unique(obs$cell_id_1km)[u]
    
    # Select data from a single sampling unit
    o = obs[obs$cell_id_1km == su,]
    
    for(x in 1:max(o$seq)){ #repeat for each sequence 
      
      # Select the sequence (i.e. column)
      indx = seq(from = 1, to = max(o$seq), by = 1)[x]
      
      # select num active cams
      n = obs$num_cams_active_at_date[obs$cell_id_1km == su & obs$seq == indx] ## COME HERE IF YOU HAVE MORE OBS COVS TO ADD! 
      
      # Add conditional to force n to be 1 if no cams detected a species, but were active in the date sequence
      if(length(n) == 0 & indx <= max(o$seq)){ n = 1 } 
      
      
      ## Fill in the obs dataframe, matching per row and column
      obs.mat[su,indx] = n
      
    } # End per sequence
  } # End per sampling unit
  
  
  ### Compress observation matrix to match sampling occasion window 
  
  ## Create a new and empty compressed matrix to fit sampling occasions
  obs2 = matrix(NA, nrow = nrow(obs.mat), ncol = round(dur/w))
  
  
  for(u in 1:nrow(obs.mat)){ #Repeat for each row in the matrix
    
    for(p in 1:round(dur/w)){ # Repeat for each sampling occasion window 
      
      # Outline the start dates of sampling occasion, using values provided in for-loop
      starts<-seq(from=0, to=dur, by=w)
      
      # Select a single start date
      l<-starts[p]
      
      # Make if-else statement 
      ifelse(all(obs.mat[u,c(l:(l+w))]== "NA", na.rm=FALSE) == "TRUE", 
             #if all values in matrix @ row u across the sampling occasion window are all NA,
             obs2[u,p]<-NA, # then leave the sampling occasion as NA,
             obs2[u,p]<- sum(as.numeric(obs.mat[u,c(l:(l+w))]), na.rm = TRUE)) 
      # But if FALSE, take the the sum of observation covariate
      
    } # End loop per sampling occasion
    
  } # End loop per row in matrix 
  
  ## Standardize the observation covaraite
  scaled.obs = scale(obs2)
  
  ## Convert observation covarites into a list and assign variable name
  obs3 = list(num_cams_active_at_date = as.data.frame(scaled.obs))
  
  # Verify the meta matches the order as the Detection history matrix 
  m = m[order(match(m$cell_id_1km, rownames(dh.mat))),]
  
  # make the umf
  if(type == "abundance"){
    umf = unmarkedFramePCount(y = dh.mat, siteCovs = m, obsCovs = obs3)
  }
  
  if(type == "occu"){
    umf = unmarkedFrameOccu(y = dh.mat, siteCovs = m, obsCovs = obs3)
  }
  
  # Save it! 
  umf.list[[t]]= umf
  names(umf.list)[t]= sp
  
  # Which species was just completed?
  print(sp)
  
}
