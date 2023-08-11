# File: 20210801_confounders.R
# Date: 2021.08.01

library(dagitty)
dag1=dagitty('dag{
                  Y4[pos="4,0"]
                  Y3[pos="3,0"]
                  Y2[pos="2,0"]
                  Y1[pos="1,0"]
                  X4[pos="3.8,1"]
                  X3[pos="2.8,1"]
                  X2[pos="1.8,1"]
                  X1[pos="0.8,1"]
                  C4[pos="3.5,-0.5"]
                  C3[pos="2.5,-0.5"]
                  C2[pos="1.5,-0.5"]
                  C1[pos="0.5,-0.5"]
                  
                  X1->X2->X3->X4
                  Y1->Y2->Y3->Y4
                  X4->Y4
                  X3->Y3
                  X2->Y2
                  X1->Y1
                  X3->Y4
                  X2->Y3
                  X1->Y2
                  Y3->X4
                  Y2->X3
                  Y1->X2
                  
                  C1->C2->C3->C4
                  C1->Y1
                  C1->X1
                  C2->Y2
                  C2->X2
                  C3->Y3
                  C3->X3
                  C4->Y4
                  C4->X4
                  
                  Y1->C2
                  Y2->C3
                  Y3->C4
                  X1->C2
                  X2->C3
                  X3->C4
             }')
plot(dag1)

paths(dag1,"X4","Y4")
adjustmentSets(dag1,exposure="X4",outcome="Y4")
adjustmentSets(dag1,exposure="X3",outcome="Y4")
adjustmentSets(dag1,exposure="X2",outcome="Y4")

