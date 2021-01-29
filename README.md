# PATMO_sulfur_cycle
1. Changes made in the code added: PATMO budget module, patmo_budget.f90.\
2. In patmo_ode.f90 added yield for COS from CS2 + DMS reactions towards COS (I did however remove these reactions + yields to see if steady state S was reached)\
3. We calculated much higher emissions and added that as input \
4. Due to high emisssions we tuned the dry deposition. \
5. Added end_of_run parameter to have output only at the end of the run. \

