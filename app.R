library(shiny)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(mrgsolve)
library(minqa)
library(MASS)
library(phonTools)
library(gridExtra)
library(deSolve)
library(cowplot)
library(plotly)
library(gapminder)
library(shinyTime)
source('helpers_temp_mtx_v8.3.R')

qweek <- sprintf('Q%sWK', 1:8)
h1 <- tags$div(textInput('targetTrough', 'Target Trough (mcg/mL)', width = '150px'), style = "display:inline-block")
h2 <- tags$div(actionButton('findTarget', 'Find'), style = "display:inline-block")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Therapeutic Drug Monitoring"),
  tabsetPanel(
    tabPanel("Application",
      fluidRow(
        column(3,
          hr(),
          fileInput('file1', 'Upload Dosing Profile (CSV)'),
          hr(),
          helpText(strong("Drug: "), "High-Dose Methotrexate"),
          hr(),
          helpText(strong("Infusion Unit: "), "mg"),
          tabsetPanel(type="pills",
          tabPanel("Baseline Information",
                   textInput("bsa", "BSA (m²)", value = "1.97"),
            selectInput("sex", "Sex:",
                        c("Male", "Female")),
            fluidRow(
              column(12, textInput("SCR_mgdl1", "Baseline Serum Creatinine (mg/dL)", value = "0.77")) ,
              column(12,shinyjs::hidden(airDatepickerInput(
                "SCR_time1",
                label = "Time",
                timepicker = TRUE,
                dateFormat = "yyyy-MM-dd",
                value = Sys.time() ) ))
            ),
          ),          
          tabPanel("Dosing",
            htmlOutput('inf1')
          ),
          tabPanel("Drug Levels",
                   htmlOutput('conc')
          ),
          tabPanel("Serum Creatinine Levels Post-infusion",
                   fluidRow(
                     column(12, textInput("SCR_mgdl2", "Serum Creatinine 2 (mg/dL)", value = "0.80")),
                     column(12,airDatepickerInput(
                       "SCR_time2",
                       label = "Time",
                       timepicker = TRUE,
                       dateFormat = "yyyy-MM-dd",
                       value = Sys.time() + 86400,
                       width = '180px'))
                   ),                   fluidRow(
                     column(12, textInput("SCR_mgdl3", "Serum Creatinine 3 (mg/dL)", value = "")),
                     column(12,airDatepickerInput(
                       "SCR_time3",
                       label = "Time",
                       timepicker = TRUE,
                       dateFormat = "yyyy-MM-dd",
                       value = NULL) )
                   ),                   fluidRow(
                     column(12, textInput("SCR_mgdl4", "Serum Creatinine 4 (mg/dL)", value = "")),
                     column(12,airDatepickerInput(
                       "SCR_time4",
                       label = "Time",
                       timepicker = TRUE,
                       dateFormat = "yyyy-MM-dd",
                       value = NULL) )
                   ),                   fluidRow(
                     column(12, textInput("SCR_mgdl5", "Serum Creatinine 5 (mg/dL)", value = "")),
                     column(12,airDatepickerInput(
                       "SCR_time5",
                       label = "Time",
                       timepicker = TRUE,
                       dateFormat = "yyyy-MM-dd",
                       value = NULL) )
                   ))
          )
        ),
        column(6,
          tags$br(),
          actionButton('runmodel', 'Create Output'),
          plotlyOutput("responseplot", width = "60vw", height = "70vh"),
          align = 'center',
          helpText(strong("Recommendation:")), 
          helpText("Follow population prediction line until first observed drug level \n then follow individual MTXPK.org prediction line"),
          hr(),
        )
      )
    ),
    tabPanel("Documentation",
             h1("Functionality"),
             p("This application takes as input patient characteristics/dosing schedule and presents the population-level expected response 
               as well as the individual-level expected response. The individual-level responses are rendered with empirical Bayesian estimates (EBES),
               so they require an observed blood level draw. The Bayesian process is based on the minimization of the likelihood with respect to the individual
               random effects [1]. The tool uses two models: Taylor et al. [2] and Blackman et al. [3]. This is based on the findings in Blackman et al., a validation
               study of high-dose methotrexate (HD-MTX) models. The study found that Blackman et al. was superior at high, early concentrations where toxicities are most likely,
               and that Taylor et al. was superior at low, later concentrations (especially with the use of EBEs. The unique ability of this application to provide
               predictions from both models enables the clinician to make the best decision at any time during a patient's concentration-time curve."),
             h1("Resources"),
             HTML("<p>The models implemented in this application are all available in the scientific literature.
                  The calculation of the EBEs for individual random effects is outlined in
                  <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/'>Kang 2012</a> [1]
                  and the code to calculate the EBEs is adapted from the open source TDM software 
                  <a href='https://mrgsolve.org/'>mrgsolve</a>.
                  </p>"),
             h1("Bibliography"),
             tags$ol(
               tags$li("Kang, D., Bae, K.-S., Houk, B. E., Savic, R. M. & Karlsson, M. O. Standard Error of Empirical Bayes Estimate in NONMEM® VI. Korean J Physiol Pharmacol (2012)."),
               tags$li("Taylor et al. MTXPK.org: A Clinical Decision Support Tool Evaluating High-Dose Methotrexate Pharmacokinetics to Inform Post-Infusion Care and Use of Glucarpidase, Clinical Pharmacology and Pharmacometrics (2020)."),
               tags$li("Blackman et al. Development and Validation of High-Dose Methotrexate Population Pharmacokinetic Models to Inform Clinical Decisions on Dosing, Journal Name (2026)."))
    ),
    tabPanel("View Profile",
      DT::dataTableOutput('profile1'),
      downloadButton('downloadTable1', 'Download Data (CSV)')
    ),
    tabPanel("Predictions",
    DT::dataTableOutput('profile2'),
    downloadButton('downloadTable2', 'Download Data (CSV)')
    )
  )
)

data2Inputs <- function(dd, myid) {
  doseLabel <- c('Bolus','Infusion','Conc')[match(substr(myid, 1, 1), c('b','i','c'))]
  nn <- nrow(dd)
  myclick <- sprintf('Shiny.onInputChange( \"delete_%s_button\" , this.id, {priority: \"event\"})', myid)
  divstyle <- "display:inline-block;vertical-align:top;height:40px"
  inputs <- vector('list', nn + 2)
#   ids <- paste0('delete_', dd[,'lab'])

  for(i in seq_len(nn)) {
    dose <- textInput(paste0('dose_', myid, '_', i), label = NULL, value = as.character(dd[i,'dose']), width = '60px')
    duration <- textInput(paste0('duration_', myid, '_', i), label = NULL, value = as.character(dd[i,'duration']), width = '60px')
    if(doseLabel == 'Conc') {
      divduration = NULL
    }else{
      divduration = tags$div(duration, style = divstyle)
    }
    suppressMessages(
      dt <- shinyWidgets::airDatepickerInput(paste0('dt_', myid, '_', i), label = NULL, value = dd[i,'dt'], timepickerOpts = list(timeFormat = 'HH:mm'),
                                             timepicker = TRUE, width = '180px', update_on = 'close'
      )
    )
#     ab <- actionButton(ids[i], label = NULL, icon = icon('trash'), onclick = myclick)
    inputs[[i + 1]] <- tags$div(tags$div(dose, style = divstyle), tags$div(dt, style = divstyle), divduration)
  }
  adder <- sprintf('Shiny.onInputChange( \"add_%s_button\" , this.id, {priority: \"event\"})', myid)
  add_label <- 'add dose'
  if(doseLabel == 'Conc') add_label <- 'add conc'
  addit <- actionButton('addit', label = add_label, icon = icon('plus'), onclick = adder, style = 'padding:4px;font-size:80%')
  h1 <- tags$div(doseLabel, style = "display:inline-block;vertical-align:top;width:60px")
  h2 <- tags$div('Time', style = "display:inline-block;vertical-align:top;width:180px")
  h3 <- tags$div('Duration', style = "display:inline-block;vertical-align:top;width:60px")
  
  #   h3 <- tags$div(actionButton('addit', label = 'add', icon = icon('plus'), onclick = adder, style = 'padding:2px;font-size:80%'), style = "display:inline-block;vertical-align:top")
  if(doseLabel == 'Conc'){h3 = NULL}
  inputs[[1]] <- tags$div(h1, h2, h3)
  inputs[[nn + 2]] <- tags$div(addit)
  inputs
}

server <- function(input, output) {
  v <- reactiveValues(dat = NULL, plot1 = NULL, plot2 = NULL, plot3 = NULL, params = NULL, sched = NULL)
  
  formulaText <- reactive({
    paste(input$drug)
  })

  # Return the formula text for printing as a caption ----
  output$caption <- renderText({
    formulaText()
  })

  observeEvent(input$file1, {
    if(!is.null(input$file1)) {
      inFile <- input$file1
      if(!is.null(inFile) && inFile$type %in% c('text/csv', 'text/comma-separated-values', 'text/plain')) {
        v$dat <- read.csv(inFile$datapath)
      } else {
        v$dat <- NULL
      }
    }
  })
  PKprof <- reactiveValues(p = NULL, Omega = NULL, wt = NULL, unit = 'mg')
  observe({
    toggleElement(id = "wt", condition = is.null(v$dat$Weight))
    toggleElement(id = "downloadTable", condition = is.null(input$file1))
    toggleElement(id = "downloadPlot", condition = !is.null(v$params))
    toggleElement(id = "vp", condition = !is.null(v$params))
  })

  ####################    ####################    ####################    ####################    ####################    ####################
 # def_time1 <- as.POSIXct(format(Sys.time(), "%Y-%m-%d"))
  def_time1 <- Sys.time() #-7200
  
  def_inf_seq <- seq(def_time1, length.out=20, by= 1800 ) # HERE 3600*24*7*2: every 2 weeks ; 3600*12: 12 hour increments
  def_inf_seq1 <- def_inf_seq[1] # HERE 3600*24*7*2: every 2 weeks ; 3600*12: 12 hour increments


  values <- reactiveValues(
    infusion1dat = data.frame(dose = c(2330), dt = def_inf_seq1, lab = 1, duration = c(24), valid = TRUE), # HERE change default
    concdat = data.frame(dose = 20, dt = def_time1 + 86400, lab = 1, valid = TRUE), # HERE; 1 hr before next dose for obs. conc time
    init = FALSE
  )
  
  
  factory_add <- function(dat, myid) {
    button <- sprintf('add_%s_button', myid)
    observeEvent(input[[button]], {
      if(!values$init) return()
      
      ix <- which(values[[dat]][,'valid'])
      
      if(length(ix) == 0) {
        # If no doses yet, start at current hour
        next_time <- as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:00"))
        first_dur <- 0
      } else {
        # Use first dose's time + first dose's duration (in hours) for next_time
        first_time <- as.POSIXct(values[[dat]][1, 'dt'])
        first_dur <- as.numeric(values[[dat]][1, 'duration'])
        next_time <- first_time + first_dur * 3600  # duration in seconds
      }
      
      # Compute duration for new row: 24 minus first dose duration
      new_duration <- 24 - first_dur
      if(new_duration < 0) new_duration <- 0  # prevent negative
      
      # Add new row
      values[[dat]] <- rbind(
        values[[dat]],
        data.frame(
          dose = NA,
          dt = next_time,
          lab = nrow(values[[dat]]) + 1,
          duration = new_duration,
          valid = TRUE
        )
      )
    })
  }
  
  
  
  factory_addconc <- function(dat, myid) {
    button <- sprintf('add_%s_button', myid)
    observeEvent(input[[button]], {
      if(!values$init) return()
      ix <- which(values[[dat]][,'valid'])
      if(is.na(ix[1])) {
        next_time <- as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:00"))
      } else {
        next_time <- as.POSIXct(format(max(values[[dat]][ix,'dt']), "%Y-%m-%d %H:00")) + 3600 * 12
      }
      values[[dat]] <- rbind(values[[dat]], data.frame(dose = NA, dt = next_time, lab = nrow(values[[dat]]) + 1, valid = TRUE))
    })
  }
  factory_del <- function(dat, myid) {
#     button <- sprintf('delete_%s_button', myid)
#     observeEvent(input[[button]], {
#       if(!values$init) return()
#       selectedId <- as.numeric(strsplit(input[[button]], "_")[[1]][2])
#       selectedRow <- match(selectedId, values[[dat]][,'lab'])
#       values[[dat]][selectedRow, 'valid'] <- FALSE
#     })
  }
  factory_dose <- function(dat, myid) {
    doseid <- sprintf('dose_%s_', myid)
    observeEvent({jj <- nrow(values[[dat]]); lapply(paste0(doseid, seq(jj)), function(i) input[[i]])}, {
      if(!values$init) return()
      labix <- values[[dat]][,'lab']
      curr_dose <- values[[dat]][,'dose']
      dose_vals <- character(length(labix))
      for(i in seq_along(dose_vals)) {
        tmp <- input[[paste0(doseid, i)]]
        if(is.null(tmp)) tmp <- NA
        dose_vals[i] <- tmp
      }
      dose_ins_row <- match(seq_along(labix), labix)
      curr_dose[dose_ins_row] <- as.numeric(dose_vals)
      values[[dat]][,'dose'] <- curr_dose
    })
  }
  factory_time <- function(dat, myid) {
    timeid <- sprintf('dt_%s_', myid)
    observeEvent({jj <- nrow(values[[dat]]); lapply(paste0(timeid, seq(jj)), function(i) input[[i]])}, {
      if(!values$init) return()
      labix <- values[[dat]][,'lab']
      curr_time <- unclass(values[[dat]][,'dt'])
      time_vals <- vector('list', length(labix))
      for(i in seq_along(time_vals)) {
        tmp <- input[[paste0(timeid, i)]]
        if(is.null(tmp)) tmp <- NA
        time_vals[[i]] <- tmp
      }
      time_ins_row <- match(seq_along(labix), labix)
      curr_time[time_ins_row] <- unlist(time_vals)
      curr_time <- as.POSIXct(curr_time, origin = '1970-01-01 00:00:00')
      values[[dat]][,'dt'] <- curr_time
    })
  }
  factory_ui <- function(dat, myid) {
    button1 <- sprintf('add_%s_button', myid)
#     button2 <- sprintf('delete_%s_button', myid)
    renderUI({
      if(!values$init) values$init <- TRUE
      ix <- isolate(c(TRUE, values[[dat]][,'valid'], TRUE))
      xx <- isolate(data2Inputs(values[[dat]], myid))
      # only add/delete should re-render
      input[[button1]]
#       input[[button2]]
      if(!is.null(xx)) do.call(tags$div, xx[ix])
    })
  }
  factory_duration <- function(dat, myid) {
    durationid <- sprintf('duration_%s_', myid)
    observeEvent({jj <- nrow(values[[dat]]); lapply(paste0(durationid, seq(jj)), function(i) input[[i]])}, {
      if(!values$init) return()
      labix <- values[[dat]][,'lab']
      curr_dur <- values[[dat]][,'duration']
      dur_vals <- character(length(labix))
      for(i in seq_along(dur_vals)) {
        tmp <- input[[paste0(durationid, i)]]
        if(is.null(tmp)) tmp <- NA
        dur_vals[i] <- tmp
      }
      dur_ins_row <- match(seq_along(labix), labix)
      curr_dur[dur_ins_row] <- as.numeric(dur_vals)
      values[[dat]][,'duration'] <- curr_dur
    })
  }
  

  
  
  factory_add('infusion1dat', 'i1')
  factory_del('infusion1dat', 'i1')
  factory_dose('infusion1dat', 'i1')
  factory_time('infusion1dat', 'i1')
  factory_duration('infusion1dat', 'i1')
  
  factory_addconc('concdat', 'c')
  factory_del('concdat', 'c')
  factory_dose('concdat', 'c')
  factory_time('concdat', 'c')

  output$inf1 <- factory_ui('infusion1dat', 'i1')
  output$conc <- factory_ui('concdat', 'c')
  ####################    ####################    ####################    ####################    ####################

  
# input and factory and data2inputs are areas that need change
# hardcoded vectors of levels + times , then need java capability to add rows
# factory was originally for multiple dosing
# start w simple case of y amount of rows of scr 
  
#   observe({
  observeEvent(input$runmodel, {
    def_start <- as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:%M"))
    def_inf <- data.frame(dose = 2330, duration = 24, dt = def_start)
    def_con <- data.frame(dose = 20, dt = def_start + 86400) # 4.5 hours
    inf_duration <- 24
    in_infusmat1 <- checkInputMatrix(values$infusion1dat[values$infusion1dat[,'valid'],c("dose","duration", "dt")], def_inf, inf_duration)
    in_concmat <- checkInputMatrix(values$concdat[values$concdat[,'valid'],c("dose","dt")], def_con)

    #repair broken inputs
    bsa <- as.numeric(input$bsa)
    #SCR_mgdl <- as.numeric(input$SCR_mgdl2)
    pt_gender <- input$sex
    negNA <- function(x) is.na(x) || x < 0
    if(negNA(bsa)) {
      bsa <- 1.97
      updateTextInput(getDefaultReactiveDomain(), inputId = "bsa", value = bsa)
    }
    scr_vals <- sapply(1:5, function(i) {
      x <- input[[paste0("SCR_mgdl", i)]]
      if (is.null(x) || x == "") NA_real_ else as.numeric(x)
    })
    

    scr_times <- sapply(1:5, function(i) {
      t <- input[[paste0("SCR_time", i)]]
      if (is.null(t)) NA else t })
    

    if(negNA(scr_vals[1])) {
      scr_vals[1] <- 68.08/88.4
      updateTextInput(getDefaultReactiveDomain(), inputId = "SCR_mgdl1", value = scr_vals[1] )
    }
    
    if(negNA(pt_gender)) {
      pt_gender <- "Male"
      updateTextInput(getDefaultReactiveDomain(), inputId = "pt_gender", value = pt_gender)
    }
    
    myparams <- list(in_infusmat1 = in_infusmat1, in_concmat = in_concmat, 
                               drug=input$drug, bsa=bsa, scr_vals=scr_vals,  scr_times = scr_times, pt_gender=pt_gender)
    # print(myparams$in_infusmat1)
    # print(myparams$in_concmat)
    # print(scr_vals)
    # print(scr_times)
    
    v$params <- myparams
    
    #set up model from helpers
    current_stuff <- reactive({
      do.call(setupModel, myparams)
    })
    
   stuff <- current_stuff() 
    v$sched <- list(stuff[[1]], stuff[[2]])
    PKprof$Omega <- stuff[[3]]$Omega
    PKprof$p <- stuff[[3]]
    #PKprof$wt <- stuff[[5]]$wt
    
   output$responseplot <- renderPlotly({
     df = stuff[[2]]
     gluc = data.frame(time = c(24, 36, 42, 48), amt = c(50, 30, 10, 5)*0.454 ) 
     concdat = stuff[[1]]
     p <- makePlots(stuff[[2]], stuff[[1]]) 
     gp <- ggplotly(p)
     con1 = df[,'Conc'] ; con2 = gluc$amt ; con3 = concdat[,'y'] ; con4 = df[,'SCR_mmol']/88.4
     myvec = c(unlist(con1[,1]), con2, unlist(con3), unlist(con4[,1]))
     y_range <- c(-0.5, max(myvec) + 2)
     y_buffer <- (y_range[2] - y_range[1]) * 0.05
     
     ceiling = 5 * round(y_range[2] / 5)
     tick_positions <- if (ceiling >= 30) {
       seq(0, ceiling, by = 10)} 
     else {
       seq(0, ceiling, by = 5)
     }
     tick_labels <- round(tick_positions / 30, 2)
     plot_range <- c(-0.5, y_range[2] * 1.05)
     
     for(i in seq_along(gp$x$data)) {
       if(grepl("Serum Creatinine", gp$x$data[[i]]$name)) {
         gp$x$data[[i]]$yaxis <- "y2"
       }
     }     
     # 1. CLEAN NAMES, MANAGE VISIBILITY, AND TARGET LINE WIDTHS
     for(i in seq_along(gp$x$data)) {
       # Clean up the name string
       raw_name <- gp$x$data[[i]]$name
       clean_name <- gsub("^\\(|\\,1\\)$", "", raw_name)
       
       if (clean_name == " " || is.null(clean_name)) {
         gp$x$data[[i]]$name <- " "
         gp$x$data[[i]]$showlegend <- TRUE   # Force it to show!
         gp$x$data[[i]]$visible <- TRUE      # Ensure it's not 'legendonly'
         gp$x$data[[i]]$marker$opacity <- 0  # Make the icon invisible
         gp$x$data[[i]]$hoverinfo <- "none"  # Don't show hover for the blank
         next 
       }
       
       # FIX: If the trace name is empty, "Trace 9", or "df_ribbon", hide it from legend
       if (is.null(clean_name) || clean_name == "" || grepl("Trace", clean_name) || grepl("df_ribbon", clean_name)) {
         gp$x$data[[i]]$showlegend <- FALSE
       } else {
         gp$x$data[[i]]$name <- clean_name
         gp$x$data[[i]]$showlegend <- TRUE
         
         # Determine if Point or Line
         is_point_data <- clean_name %in% c("Observed Drug Level", 
                                            "Serum Creatinine (mg/L)", 
                                            "Glucarpidase Consensus Guidelines")
         
         if(is_point_data) {
           gp$x$data[[i]]$mode <- "markers"
           gp$x$data[[i]]$line$width <- 0 
         } else {
           gp$x$data[[i]]$mode <- "lines"
           # Thicken Blackman Pop/Indiv lines
           if(grepl("Blackman et al.", clean_name)) {
             if(!grepl("probability", clean_name)) {
               gp$x$data[[i]]$line$width <- 2.5    
             } else {
               gp$x$data[[i]]$line$width <- 1    
             }
           } else {
             gp$x$data[[i]]$line$width <- 1.5     
           }
         }
         
         # Handle Secondary Axis for SCr
         if(grepl("Serum Creatinine", clean_name)) {
           gp$x$data[[i]]$yaxis <- "y2"
         }
         
         gp$x$data[[i]]$legendgroup <- clean_name
       }
       
       
     }
     
     # 2. FORCE ORDER (Plotly shows traces in the order they appear in gp$x$data)
     level_order <- c(
       "Observed Drug Level",              # Col 1, Row 1
       " ",                                # Col 2, Row 1 (The Blank)
       "Serum Creatinine (mg/L)",          # Col 1, Row 2
       "Glucarpidase Consensus Guidelines", # Col 2, Row 2
       "Population level: Blackman et al.",
       "Individual level: Blackman et al.",
       "Population level: MTXPK.org",
       "Individual level: MTXPK.org",
       "Upper 95% population probability: Blackman et al.",
       "Lower 95% population probability: Blackman et al."
     )
     
     # This re-sorts the traces based on your preferred order
     trace_names <- sapply(gp$x$data, function(x) x$name)
     order_indices <- match(level_order, trace_names)
     order_indices <- order_indices[!is.na(order_indices)]
     
     # Append any traces not in our list (like the ribbon) to the end
     other_indices <- setdiff(seq_along(gp$x$data), order_indices)
     gp$x$data <- gp$x$data[c(order_indices, other_indices)]
     
    
   
     gp %>% 
       layout(
         legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2, traceorder = "normal"),
         margin = list(r=80),
         yaxis = list(range = plot_range, 
                      tickmode = "auto", 
                      tick_vals = tick_positions,
                      ticks = "outside",           # Draws the actual tick lines
                      ticklen = 5,                 # Length of the tick marks in pixels
                      tickwidth = 1,               # Thickness of the tick marks
                      tickcolor = "black"      # Match your SCr point/font color
         ),
         yaxis2 = list(overlaying = "y", 
                       side = "right", 
                       range = plot_range,
                       tickvals = tick_positions, 
                       ticktext = round(tick_positions / 30, 2),
                       title = list(text ="Serum Creatinine (mg/L)", font = list(size = 14), color = "purple"), 
                       ticktext = tick_positions,tickfont = list(size = 12, color = "purple"),
                       ticks = "outside",           # Draws the actual tick lines
                       ticklen = 5,                 # Length of the tick marks in pixels
                       tickwidth = 1,               # Thickness of the tick marks
                       tickcolor = "black"      # Match your SCr point/font color
         )
       )
   })
     
     
   #    gp %>% 
   #     layout(
   #     legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2),
   #     margin = list(r=80),
   #     yaxis = list(range = plot_range),
   #     yaxis2 = list(
   #       overlaying = "y",
   #       side = "right",
   #       range = plot_range,
   #       tickvals = tick_positions,
   #       ticktext = round(tick_positions / 30, 2),
   #       tickfont = list(size = 11, color = "grey30"),
   #       showgrid = FALSE, visible = TRUE, tickcolor="black",
   #       title = list(text ="Serum Creatinine (mg/L)", 
   #                                      font = list(size = 16, family = "Arial", color = "black"))
   #     )#,
   #    # annotations = list(text = "fngjgfk")
   #   ) #%>% 
   #   #    style( mode = "lines+markers", 
   #   #   marker = list(symbol = "diamond-open", size = 10, color = "#0000ee", width = 7),
   #   #   line = list(width = 0, color = "#0000ee"),
   #   #   traces = 8
   #   # ) 
   #   
   # })
  })





  tableInput <- function() {
    if(is.null(v$params)) return(NULL)
    
    i <- v$params[["in_infusmat1"]]
    c <- v$params[["in_concmat"]]
    w <- 1 # multiply by "weight" or not?
    c.s <- data.frame(Time = as.numeric(c[,'time'])/3600, DV = c[,1], AMT = NA, MDV = 0, EVID = 0, Rate = NA, Duration = NA, SCR = NA)
    myrate <- round(as.numeric(i[,1]) / as.numeric(i[,"duration"]),2)
    i.s <- data.frame(Time = as.numeric(i[,'time'])/3600, DV = NA, AMT = i[,1], 
                      MDV = 1, EVID = 1, Rate = myrate, Duration = i[,"duration"], SCR = NA)
    t.out <- rbind(i.s,c.s)
    zero.amt <- which(!is.na(t.out[,'AMT']) & t.out[,'AMT']==0)
    if(length(zero.amt)) {
      t.out <- t.out[-zero.amt,]
    }
    t.out <- t.out[order(t.out[,'Time']),]

    scr_vals=round(v$params$scr_vals,2) ; scr_times = v$params$scr_times/3600
    t.out[which(t.out$Time == min(t.out[,'Time'])),"SCR"] = scr_vals[1]
    
    if (length(scr_times > 1)){
    scr.s = data.frame(Time = scr_times[2:length(scr_times)], 
                       DV = NA, AMT = NA, MDV = 0, EVID = 0, Rate = NA, Duration = NA, SCR = scr_vals[2:length(scr_vals)])
    }
    
    t.out1 <- rbind(t.out, scr.s)
    t.out1 <- t.out1[order(t.out1[,'Time']),]
    t.out1[,'Time'] <- round((t.out1[,'Time'] - min(t.out[,'Time'])), 2)
    
    s = which(table(t.out1$Time) > 1)
    names = as.numeric(names(s))
    
    if (length(names) > 0){
      for (i in 1:length(names)) {
        idx = which(t.out1$Time == names[i])
        DV = t.out1[idx,"DV"][!is.na(t.out1[idx,"DV"])]
        SCR = t.out1[idx,"SCR"][!is.na(t.out1[idx,"SCR"])]
        t.out1[idx[1], "DV"] = DV
        t.out1[idx[1], "SCR"] = SCR
        t.out1 = t.out1[-idx[-1],]
        }
    }
    

    t.out1$bsa <- unique(v$params$bsa)
    t.out1$Sex <- unique(v$params$pt_gender)
    t.out1 %>%  filter(!is.na(Time)) %>% rename(`BSA (m²)`= bsa, `SCR (mg/L)` = SCR) %>%  
      fill(`SCR (mg/L)`, .direction = "down") #%>%  
     # fill(`BSA (m²)`, .direction = "down")%>%  
     # fill(Sex, .direction = "down")
  }
  
  table2Input <- function() {
    if(is.null(v$sched[[2]])) return(NULL)
    
    d <- v$sched[[2]]
    #max = ceiling(max(d$time)/6)
    dat = d %>% 
      filter(time %in% seq(from = 0, to = 600, by = 6)) %>% 
      dplyr::select(time, bsa, SCR_mmol, pt_gender, Conc, tag) %>% 
      mutate(pt_gender = ifelse(pt_gender == 0, "Male", "Female"),
             Conc = round(Conc, 2),
             SCR_mmol = round(SCR_mmol/88.4, 2)) %>% 
      rename(`BSA (m²)` = bsa,
             `SCR (mg/L)` = SCR_mmol,
             `Concentration (mg/L)` = Conc,
             Sex = pt_gender,
             `Time (hrs)` = time)
    
    dat <- dat %>%
      pivot_wider(
        id_cols = c(`Time (hrs)`, `BSA (m²)`, `SCR (mg/L)`, Sex),
        names_from = tag,
        values_from = `Concentration (mg/L)`
      )
    
    dat
  }
  
  output$profile1 <- DT::renderDT({
    tableInput()
  })
  
  output$profile2 <- DT::renderDT({
    table2Input()
  })


  output$downloadPlot <- downloadHandler(
    filename = function() { sprintf('TDM_plot_%s.pdf', Sys.Date()) },
    content = function(file, width = 8, height = 6) {
      pdf(file, width = width, height = height, onefile = TRUE)
      grid.arrange(v$plot1, v$plot2, v$plot3, heights = c(5, 4, 4))
      dev.off()
  })

  output$downloadTable <- downloadHandler(
    filename = function() { sprintf('TDM_profile_%s.csv', Sys.Date()) },
    content = function(file) {
      write.csv(tableInput(), file, row.names = FALSE)
  })

  output$downloadTable1 <- downloadHandler( 
    filename = function() { sprintf('TDM_profile_%s.csv', Sys.Date()) },
    content = function(file) {
      write.csv(tableInput(), file, row.names = FALSE)
  })
  
  output$downloadTable2 <- downloadHandler( 
    filename = function() { sprintf('TDM_profile_%s.csv', Sys.Date()) },
    content = function(file) {
      write.csv(table2Input(), file, row.names = FALSE)
    })

  output$imgdat <- renderUI({
    tags$img(src = "https://raw.githubusercontent.com/michaelleewilliams/michaelleewilliams.github.io/master/pictures/datex.png",height="50%", width="50%", align="center")
  })

  output$imgscn <- renderUI({
    tags$img(src = "https://raw.githubusercontent.com/michaelleewilliams/michaelleewilliams.github.io/master/pictures/scn.png",height="100%", width="100%", align="center")
  })
}

shinyApp(ui, server)
