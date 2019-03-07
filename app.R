#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)
library(ggplot2)
library(gridExtra)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Population growth activity: harvesting"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(position='left',
      sidebarPanel(
         sliderInput("r","r, Population growth rate",
                     min = -0.5,
                     max = 1.5,
                     value = 0.8),
         sliderInput("K","K, Carrying capacity",
                     min = 1,
                     max = 100,
                     value = 50),
         sliderInput("N0","N(0), Initial abundance",
                     min = 0,
                     max = 100,
                     value = 1),
         checkboxInput("allow.fixed.H", "Allow fixed harvesting?", value = F),
         sliderInput("H","H, fixed harvest rate",
                     min = 0,
                     max = 0.25*1.5*100+0.1,
                     value = 0)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        column(6,plotOutput(outputId="popPlot", width="500px",height="400px")) 
        #plotOutput("popPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
  ### Run background calculations
  
  plot1<-reactive({
    if(!input$allow.fixed.H){
    
      # Differential equation set up
      Logistic<-function(t,state,parameters){
        with(as.list(c(state,parameters)),{
          # rate of change
          dN<- r*N*(1-N/K)
        
          # return the rate of change
          list(c(dN))
        })	# end of with(as.list...
      }
    
      # define parameters and IC's
      parameters<-c(r=input$r,K=input$K)
      state<-c(N=input$N0)
      times<-seq(0,100,0.01)
    
      # Solve ODE
      out<-ode(y=state,times=times,func=Logistic,parms=parameters,method='ode45')
      out<-data.frame(out)
    
      g1<-ggplot(out,aes(x=time,y=N))+
        geom_line(colour='blue',size=1)+
        scale_y_continuous('Abundance, N',limits=c(0,110))+
        scale_x_continuous('Time, t',limits=c(0,100))+
        theme_bw()+
        ggtitle('Population dynamics')
    }else{ # FIXED HARVEST
      
      # Differential equation set up
      LogisticH<-function(t,state,parameters){
        with(as.list(c(state,parameters)),{
          # rate of change
          dN<- r*N*(1-N/K)-H
          dY<-H
          
          # return the rate of change
          list(c(dN,dY))
        })	# end of with(as.list...
      }
      
      # define parameters and IC's
      parameters<-c(r=input$r,K=input$K,H=input$H)
      state<-c(N=input$N0,Y=0)
      times<-seq(0,100,0.01)
      
      # Solve ODE
      out<-ode(y=state,times=times,func=LogisticH,parms=parameters,method='ode45')
      out<-data.frame(out)
      
      # figure out where, if anywhere, population crashes...
      tmp1<-unlist(na.omit(out$time[out$N<0]))
      if(length(tmp1)>0){
        crash.time<-min(tmp1)
        out$Y<-ifelse(out$time<crash.time,out$Y,out$Y[out$time==crash.time])
      }

      # Total yield
      ty<-out$Y[nrow(out)]
      
      g1a<-ggplot(out,aes(x=time,y=N))+
        geom_line(colour='blue',size=1)+
        scale_y_continuous('Abundance, N')+
        scale_x_continuous('Time, t',limits=c(0,100))+
        coord_cartesian(ylim = c(0,110))+
        theme_bw()+
        ggtitle('Population dynamics')
      g1b<-ggplot(out,aes(x=time,y=Y))+
        geom_line(colour='red',size=1)+
        scale_y_continuous('Cumulative yield')+
        scale_x_continuous('Time, t',limits=c(0,100))+
        theme_bw()+
        ggtitle(paste('Final yield = ',round(ty,2)))
      g1<-grid.arrange(g1a,g1b,nrow=1)
    }
    g1
  })
    
  plot2<-reactive({
    if(!input$allow.fixed.H){
      maxK<-100
      xs<-seq(0,maxK,0.01)
      ys<-input$r*xs*(1-xs/input$K)
      g2<-ggplot(data.frame(xs,ys),aes(x=xs,y=ys))+
        geom_line(size=1,colour='blue')+
        geom_hline(yintercept = 0)+
        scale_y_continuous('dN/dt')+
        scale_x_continuous('N',limits=c(0,maxK))+
        coord_cartesian(ylim=c(0,0.25*1.5*maxK))+
        theme_bw()
    }else{#FIXED HARVEST
      maxK<-100
      xs<-seq(0,maxK,0.01)
      ys<-input$r*xs*(1-xs/input$K)
      g2<-ggplot(data.frame(xs,ys),aes(x=xs,y=ys))+
        geom_line(size=1,colour='blue')+
        geom_hline(yintercept = 0)+
        geom_hline(yintercept=input$H,colour='red',linetype=2)+
        geom_text(aes(x=c(75),y=c(input$H+5),label='Harvest rate'),colour='red')+
        scale_y_continuous('dN/dt')+
        scale_x_continuous('N',limits=c(0,maxK))+
        coord_cartesian(ylim=c(0,0.25*1.5*maxK))+
        theme_bw()
      
    }
    
    g2
  })
  
  
  output$popPlot <- renderPlot({
    plot.list<-list(plot1(),plot2())
    grid.arrange(grobs=plot.list)   
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

