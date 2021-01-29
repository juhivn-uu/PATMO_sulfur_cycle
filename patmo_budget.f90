module patmo_budget
        use patmo_commons
         implicit none
         integer,parameter::nprocess = 8
         integer,parameter::patmo_npr_chemprod = 1  
         integer,parameter::patmo_npr_chemloss = 2 
         integer,parameter::patmo_npr_diffup = 3
         integer,parameter::patmo_npr_diffdown = 4
         integer,parameter::patmo_npr_emission = 5
         integer,parameter::patmo_npr_drydep = 6
         integer,parameter::patmo_npr_wetdep = 7
         integer,parameter::patmo_npr_photolysis = 8

         integer,parameter,dimension(nprocess)::indexProcess3 = (/&
         patmo_npr_chemprod,&
         patmo_npr_chemloss,&
         patmo_npr_diffup,&
         patmo_npr_diffdown,&
         patmo_npr_emission,&
         patmo_npr_drydep,&
         patmo_npr_wetdep,&
         patmo_npr_photolysis/)

         real*8::budget(cellsNumber,speciesNumber,nprocess)
 end module patmo_budget
