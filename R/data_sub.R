################################################################################
## Purpose: Functions for substituting IHME data into EPP data object
## Date created: 
## Date modified:
## Author: Austin Carter, aucarter@uw.edu
## Run instructions: 
## Notes:
################################################################################

### Setuprm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/HIV/")

### Paths
inputs.dir <- paste0(root,"/WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/")
aim.dir <- paste0(inputs.dir, "AIM_assumptions/")

### Code
extend.trans.params <- function(dt, start.year, stop.year) {
	for(n in names(dt)) {
		## create time series of on/off art mortality and cd4 progression
		n.years <- length(start.year:stop.year)
		# mortality
		single.year.cd4artmort <- attr(dt[[n]], "eppfp")$cd4artmort
		time.series.cd4artmort <- do.call("rbind", rep(list(single.year.cd4artmort), n.years))
		attr(dt[[n]], "eppfp")$cd4artmort <- time.series.cd4artmort
		#cd4 progression
		single.year.cd4prog <- attr(dt[[n]], "eppfp")$cd4prog
		time.series.cd4prog <- do.call("rbind", rep(list(single.year.cd4prog), n.years))
		attr(dt[[n]], "eppfp")$cd4prog <- time.series.cd4prog

		attr(dt[[n]], "eppfp")$mortyears <- n.years
		attr(dt[[n]], "eppfp")$cd4years <- n.years

	}
	return(dt)
}

sub.prev.survey <- function(dt, loc) {
	supplement_survey <- fread(paste0(root,"/WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/supplement_survey_data.csv"))
	if (loc %in% supplement_survey[,unique(iso3)] & !(loc %in% c("MWI", "ZMB"))) {
	  print('found survey')
	  survey_subpop <- supplement_survey[iso3==dt_iso3, unique(subpop)]
	  print(names(dt))
	  print(survey_subpop)
	  tmp_survey <- supplement_survey[iso3==dt_iso3,.(year, prev, se, n)]
	  tmp_survey[,used:=TRUE]
	  tmp_survey[prev==0,used:=FALSE]
	  # tmp_survey[year == 2000, used:=FALSE]
	  tmp_survey[,W.hhs:=qnorm(prev)]
	  tmp_survey[,v.hhs:=2*pi*exp(W.hhs^2)*se^2]
	  tmp_survey[,sd.W.hhs := sqrt(v.hhs)]
	  tmp_survey[,idx := year - (start.year-1)]

	  attr(dt[[survey_subpop]], 'likdat')$hhslik.dat <- rbind(attr(dt[[survey_subpop]], 'likdat')$hhslik.dat, as.data.frame(tmp_survey[used==TRUE,]))
	}

}

sub.anc <- function(dt, loc) {

}

sub.pop.params <- function(dt, loc) {
	location_id <- loc.table[ihme_loc_id==ihme_loc, location_id]
    ### ZAF read in national level country file
    if ("ZAF" %in% gsub('\\_.*', "", ihme_loc)) {
    location_id <- loc.table[ihme_loc_id=="ZAF", location_id]
    } 

    gbd.pop <- get_population(year_id = -1, location_id = location_id, sex_id = 3, age_group_id = 24, gbd_round_id = 5)
    epp.adult.pop <- data.table(year_id=epp.subp$total$year, epp_pop=epp.subp$total$pop15to49)
    merged.pop <- merge(epp.adult.pop, gbd.pop[,.(year_id, population)], by="year_id")
    merged.pop[,scale:=population/epp_pop]
    scale <- merged.pop[,scale]
    scale <- c(scale, rep(scale[length(scale)], (nrow(epp.subp$total)- length(scale))))

    for (pop in names(epp.subp)) {
    temp.dt <- epp.subp[[pop]]
        if(pop=="subpops") {
            for (subpop in names(temp.dt)){
               temp.dt1 <- temp.dt[[subpop]]          
            for (col_id in c("pop15to49", "pop15", "pop50", "netmigr")){
              epp.subp[[pop]][[subpop]][,col_id] <- temp.dt1[,col_id]*scale
             }
            } 
        } else {
            for (col_id in c("pop15to49", "pop15", "pop50", "netmigr")){
                epp.subp[[pop]][,col_id] <- temp.dt[,col_id]*scale
            }
        }
    }
}

calc.expand.pop <- function(loc, sex.agg = T) {
	
	loc.id <- loc.table[ihme_loc_id==loc, location_id]
	
	## Load central functions
	if(!"get_population" %in% ls()) {
		source(paste0(root, "temp/central_comp/libraries/current/r/get_population.R"))
	}
	if(!"get_age_map" %in% ls()) {
		source(paste0(root, "Project/Mortality/shared/functions/get_age_map.r"))
		age.table <- data.table(get_age_map())
	}



	IRR2 <- fread(paste0(aim.dir, "sex_age_pattern/age_IRRs/Feb17/GEN_IRR.csv"))
	IRR2 <- IRR2[age < 55,]

	sex_IRR <- fread(paste0(aim.dir, "sex_age_pattern/FtoM_inc_ratio_epidemic_specific.csv"))
	sex_IRR <- sex_IRR[epidemic_class=="GEN",]
	sex_IRR[,year:=year+start.year-1]

	missing_years <- c()
	if (sex_IRR[,max(year)] < stop.year)
	  missing_years <- (sex_IRR[,max(year)]+1):stop.year
	replace_IRR <- sex_IRR[order(year)][rep(nrow(sex_IRR), times=length(missing_years))]
	if (length(missing_years) > 0)
	  replace_IRR[,year:=missing_years]
	sex_IRR <- rbind(sex_IRR, replace_IRR)

	sex_IRR[,sex:=2]

	male_IRR <- copy(sex_IRR)
	male_IRR[,FtoM_inc_ratio:=1.0]
	male_IRR[,sex:=1]

	sex_IRR <- rbind(sex_IRR, male_IRR)

	## Read in population for age-sex structure for aggregation
	in.pop <- get_population(age_group_id = 8:15, location_id = loc.id, year_id = -1, sex_id = 1:2, location_set_id = 21)
	pop <- merge(in.pop, age.table[, .(age_group_id, age_group_name_short)], by = "age_group_id", all.x = T)
	setnames(pop, 
	c("year_id", "sex_id", "age_group_name_short", "population"),
	c("year", "sex", "age", "value")
	)
	pop[, (setdiff(names(pop), c("year", "sex", "age", "value"))) := NULL]


	pop$age <- strtoi(pop$age)
	pop[(age-5) %%  10 != 0, age:=as.integer(age-5)]
	pop[,value:=as.numeric(value)]

	pop1 <-data.table(aggregate(value ~ sex + age + year,pop,FUN=sum))[order(sex,age)]
	missing_years <- c()
	if (pop1[,max(year)] < stop.year)
	  missing_years <- (pop1[,max(year)]+1):stop.year
	replace_pop <- pop1[rep(which(pop1[,year] == pop1[,max(year)]), times=length(missing_years))]
	replace_years <- rep(1:(stop.year-pop1[,max(year)]), each=length(which(pop1[,year] == pop1[,max(year)])))
	replace_pop[,year:=year+replace_years]
	pop1 <- rbind(pop1, replace_pop)
	IRR <- runif(16, IRR2$lower, IRR2$upper)
	IRR2[,IRR:=IRR]

	IRR2[,IRR:=IRR2[,IRR]/IRR2[age==25 & sex==1,IRR]]

	combined_IRR <- merge(sex_IRR, IRR2, by='sex', allow.cartesian=TRUE)
	combined_IRR[,comb_IRR := FtoM_inc_ratio * IRR]

	pop2 <- merge(pop1, combined_IRR, by=c('sex', 'age', 'year'))

	pop2[,wt:=comb_IRR*value]
	
	sex_agg <- pop2[,.(wt=sum(wt)),by=.(year, age)]

	total <- pop2[,.(total = sum(wt)),by=.(year)]
	pop2 <- merge(pop2, total, by=c('year'))
	pop2[,ratio:=wt/total]


	sex_agg <- merge(sex_agg, total, by=c('year'))
	sex_agg[,ratio:=wt/total]

	if(sex.agg) {
		out.pop <- sex_agg
	} else {
		out.pop <- pop2
	}

	return(out.pop)
}

sub.on.art <- function(dt, loc, k) {

}

sub.off.art <- function(dt, loc, k) {
	# Off-ART Mortality
	mortnoart <- fread(paste0(aim.dir, "transition_parameters/HIVmort_noART/current_draws/",loc,"_mortality_par_draws.csv"))
	mortnoart[,draw:=rank(-mort,ties.method="first"),by=c("age","cd4")]
	mortnoart <- mortnoart[order(age,cd4,draw)]
	mortnoart_read <- mortnoart[,c("age","cd4","draw","mort"), with=F]
	mortnoart <- mortnoart_read[draw==k,]
	mortnoart[,age:= as.integer(sapply(strsplit(mortnoart[,age],'-'), function(x) {x[1]}))]
	mortnoart[,risk:=-1*log(1-mort)/0.1]
	mortnoart[,prob:=1-exp(-1*risk)]

	cd4_cats <- unique(mortnoart[,cd4])
	cd4_vars <- data.table(cd4=cd4_cats)
	sex_agg <- calc.expand.pop(loc)
	expanded_pop <- sex_agg[rep(1:nrow(sex_agg), times=length(cd4_cats))]
	expanded_pop <- expanded_pop[order(year, age)]
	expanded_pop <- cbind(expanded_pop, cd4_vars[rep(1:(nrow(cd4_vars)), times=nrow(sex_agg)),])
	combined_mu <- merge(expanded_pop, mortnoart, by=c('cd4', 'age'))
	mortnoart <- combined_mu[,.(prob=sum(ratio*prob)), by=.(cd4, year)]

	mortnoart <- mortnoart[cd4=="GT500CD4", cat := 1]
	mortnoart <- mortnoart[cd4=="350to500CD4", cat := 2]
	mortnoart <- mortnoart[cd4=="250to349CD4", cat := 3]
	mortnoart <- mortnoart[cd4=="200to249CD4", cat := 4]
	mortnoart <- mortnoart[cd4=="100to199CD4", cat := 5]
	mortnoart <- mortnoart[cd4=="50to99CD4", cat := 6] 
	mortnoart <- mortnoart[cd4=="LT50CD4", cat := 7] 
	mortnoart[,risk:=-1*log(1-prob)]
	mortnoart <- mortnoart[,.(year,risk,cat)]
	mortnoart <- mortnoart[, setattr(as.list(risk), 'names', cat), by=.(year)]
	mortnoart <- mortnoart[order(year)]
	mortnoart <- mortnoart[,c("1","2","3","4","5","6", "7"), with=F]
	mortnoart <- data.frame(mortnoart)
	mugbd <- as.matrix(mortnoart)
	mu <- as.vector(t(mugbd))
	for (n in names(dt)) {
        attr(dt[[n]], 'eppfp')$cd4artmort[,0] <- mu
     }
	return(dt)
}


sub.cd4.prog <- function(dt, loc, k) {
	for (n in names(dt)) {
        attr(dt[[n]], 'eppfp')$cd4prog <- progdata
     }	
}

### End