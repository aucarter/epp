prep_epp_data <- function(loc, popadjust = FALSE, popupdate = FALSE, proj.end = 2017.5, stop_collapse = FALSE) {
    loc.table <- as.data.table(loc.table)
    if(stop_collapse) {
        collapse <- F
    } else {
        print(is.data.table(loc.table))
        collapse <- loc.table[ihme_loc_id == loc, collapse_subpop]
        print(is.data.table(loc.table))
        unaids.year <- loc.table[ihme_loc_id == loc, unaids_recent]
        dir <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/UNAIDS_country_data/", unaids.year, "/")
        pjnz <- paste0(dir, loc, ".PJNZ")        
    }
    if(collapse) {
        print("Collapsing subpopulations")
        epp_totals <- collapse_epp(loc)
        eppd <- epp_totals$eppd.tot
        epp.subp <- epp_totals$epp.subp.tot
        epp.input <- epp_totals$epp.input.tot
        epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)
        val <- setNames(vector("list", length(eppd)), names(eppd))
        set.list.attr <- function(obj, attrib, value.lst) mapply(function(set, value) {
            attributes(set)[[attrib]] <- value
            set
        }, obj, value.lst)
        val <- set.list.attr(val, "eppd", eppd)
        val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end))
        val <- set.list.attr(val, "country", attr(eppd, "country"))
        val <- set.list.attr(val, "region", names(eppd))

        prepped.dt <- val

    } else {
        prepped.dt <- prepare_epp_fit(pjnz)
    }
    return(prepped.dt)
}

collapse_epp <- function(loc) {
    loc.table <- as.data.table(loc.table)
    print(is.data.table(loc.table))
    unaids.year <- loc.table[ihme_loc_id == loc, unaids_recent]
    dir <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/UNAIDS_country_data/", unaids.year, "/")
    pjnz.list <- list.files(dir, pattern = "PJNZ", full.names = T)
    file.list <- grep(loc, pjnz.list, value = T)
    ## eppd
    eppd.list <- lapply(file.list, function(file) {
        pjnz <- file
        eppd <- read_epp_data(pjnz)
    }) 
    
    eppd.list <- unlist(eppd.list,recursive = FALSE )
    
    eppd.tot <- eppd.list[1]
    loc.name <- loc.table[ihme_loc_id == loc, location_name]
    subpop.tot <- paste0(loc.name," Total")
    names(eppd.tot) <- subpop.tot

    # region
    eppd.tot[[subpop.tot]]$region <- subpop.tot 
    
    #country 
    attr(eppd.tot,"country") <- eppd.tot[[1]]$country

    # anc.used (append)
    eppd.tot[[subpop.tot]]$anc.used <- unlist(lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        anc.used <- eppd$anc.used
    }))

    # anc.prev (append)
    eppd.tot[[subpop.tot]]$anc.prev <- do.call(rbind, lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        anc.prev <- eppd$anc.prev
    }))

    # anc.n (append)
    eppd.tot[[subpop.tot]]$anc.n <- do.call(rbind, lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        anc.n <- eppd$anc.n
    }))

    # hhs (append) ** be careful "not used TRUE"
    hhs.temp <- data.table(do.call(rbind, lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        hhs <- eppd$hhs
    })))
    hhs.temp <- hhs.temp[used == TRUE]
    hhs.temp[, pos := n * prev]
    hhs.sum <- hhs.temp[, lapply(.SD, sum), by = .(year)]
    hhs.sum[, prev := pos / n]
    hhs.sum[, se := ((prev * (1 - prev)) / n)**0.5]
    hhs.sum[, used := NULL]
    hhs.sum[, pos := NULL]
    hhs.sum[, used := TRUE]
    eppd.tot[[subpop.tot]]$hhs <- as.data.frame(hhs.sum)
    # eppd.tot[[subpop.tot]]$hhs <- as.data.frame(hhs.temp)

    ## epp.subp
    epp.subp.list <- lapply(file.list, function(file) {
        pjnz <- file
        epp.subp <- read_epp_subpops(pjnz)
    })
    
    #this depends on first one haveing no subpops so I think better to make an empty list
    # epp.subp.tot <- epp.subp.list[[1]]
    # names(epp.subp.tot$subpops) <- subpop.tot
    epp.subp.tot <- list()
    

    # total
    total.temp <- data.table(do.call(rbind, lapply(epp.subp.list, function(epp.subp) {
        anc.n <- epp.subp$total
    })))
    total.sum <- total.temp[, lapply(.SD, sum), by = .(year)]
    epp.subp.tot$total <- as.data.frame(total.sum)
    epp.subp.tot$subpops[[subpop.tot]] <- as.data.frame(total.sum)

    ## epp.input
    epp.input.list <- lapply(file.list, function(file) {
        pjnz <- file
        epp.subp <- read_epp_input(pjnz)
    })
    epp.input.tot <- epp.input.list[[1]]
    attr(epp.input.tot,"country") <- subpop.tot

    # start.year (check for difference)
    start.years <- unlist(lapply(epp.input.list, function(epp.input) {
        start.year <- epp.input$start.year
    }))
    length(unique(start.years)) == 1

    # stop.year (check for difference)
    stop.years <- unlist(lapply(epp.input.list, function(epp.input) {
        stop.year <- epp.input$stop.year
    }))
    length(unique(stop.years)) == 1

    # epidemic.start (check for difference)
    epidemic.starts <- unlist(lapply(epp.input.list, function(epp.input) {
        epidemic.start <- epp.input$epidemic.start
    }))
    length(unique(epidemic.starts)) == 1

    # epp.pop (sum and mean)
    epp.pop.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        epp.pop <- epp.input$epp.pop
    })))
    epp.pop.sum <- epp.pop.temp[, lapply(.SD, sum), by = .(year)]
    epp.pop.mean <- epp.pop.temp[, lapply(.SD, mean), by = .(year)]
    epp.pop.comb <- cbind(epp.pop.sum[, .(year, pop15to49, pop15, pop50, netmigr)], epp.pop.mean[, .(cd4median, hivp15yr)])
    epp.input.tot$epp.pop <- epp.pop.comb

    # cd4lowlim (check for difference)
    cd4lowlim.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        cd4lowlim <- epp.input$cd4lowlim
    })))

    # cd4initperc (check for difference)
    cd4initperc.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        cd4initperc <- epp.input$cd4initperc
    })))

    # cd4stage.dur (check for difference)
    cd4stage.dur.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        cd4stage.dur <- epp.input$cd4stage.dur
    })))

    # cd4mort, artmort.less6mos, artmort.6to12mos, artmort.after1yr (leave the same)

    # infectreduc (check for difference)
    infectreducs <- unlist(lapply(epp.input.list, function(epp.input) {
        infectreduc <- epp.input$infectreduc
    }))
    length(unique(infectreducs)) == 1

    # epp.art (sum and mean) ** beware of percentages!!! also not sure whether 1stto2ndline is count or percent
    epp.art.temp <- data.table(rbind.fill(lapply(epp.input.list, function(epp.input) {
        epp.art <- epp.input$epp.art
    })))
    epp.art.temp[is.na(m.isperc), m.isperc := "N"]
    epp.art.temp[is.na(f.isperc), f.isperc := "N"]
    if("P" %in% unique(c(epp.art.temp$m.isperc, epp.art.temp$f.isperc))) {
        pop <- epp.pop.temp[year %in% unique(epp.art.temp$year)]
        epp.art.temp <- cbind(epp.art.temp, pop[, .(pop15to49)])
        epp.art.temp[m.isperc == "P", c("m.val", "f.val") := (((m.val + f.val) / 2) / 100) * (pop15to49 / 2)]
        epp.art.temp[, pop15to49 := NULL]
    }
    epp.art.temp[, m.isperc := NULL]
    epp.art.temp[, f.isperc := NULL]
    epp.art.sum <- epp.art.temp[, lapply(.SD, sum), by = .(year)]
    epp.art.mean <- epp.art.temp[, lapply(.SD, mean), by = .(year)]
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    epp.art.mode <- epp.art.temp[, lapply(.SD, Mode), by = .(year)]
    epp.art.comb <- cbind(epp.art.sum[, .(year, m.val, f.val, artdropout)], 
                            epp.art.mode[, .(cd4thresh)],
                            epp.art.mean[, c("m.perc50plus", "f.perc50plus", "perc50plus", "1stto2ndline", "art15yr"), with = F])
    epp.art.comb[, c("m.isperc", "f.isperc") := "N"]
    epp.art.order <- epp.art.comb[, c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh", "m.perc50plus", "f.perc50plus", "perc50plus", "1stto2ndline", "art15yr"), with = F]
    epp.input.tot$epp.art <- epp.art.order

    # art.specpop (check for difference)
    art.specpop.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        art.specpop <- epp.input$art.specpop
    })))

    # hivp15yr.cd4dist (check for difference)
    hivp15yr.cd4dist.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        hivp15yr.cd4dist <- epp.input$hivp15yr.cd4dist
    })))

    # art15yr.cd4dist (check for difference)
    art15yr.cd4dist.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
        art15yr.cd4dist <- epp.input$art15yr.cd4dist
    })))

    # epidemic.type (check for difference)
    # epidemic.types <- unlist(lapply(epp.input.list, function(epp.input) {
    #   epidemic.type <- epp.input$epidemic.type
    # }))
    # length(unique(epidemic.types)) == 1

    # ## Save
    # dir.create(paste0(dir, loc), showWarnings = F)
    # save(eppd.tot, file = paste0(dir, loc, "/eppd.Rdata"))
    # save(epp.subp.tot, file = paste0(dir, loc, "/epp_subp.Rdata"))
    # save(epp.input.tot, file = paste0(dir, loc, "/epp_input.Rdata"))
    epp_totals <- list(eppd.tot = eppd.tot, epp.subp.tot = epp.subp.tot, epp.input.tot = epp.input.tot )
    return(epp_totals)

}
