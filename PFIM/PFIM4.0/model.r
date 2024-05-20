source(file.path('/Users/helingtong/Dropbox/Mac/Desktop/Tong/PFIM4.0/Program/LibraryPK.r'))

formA<-bolus_1cpt_VCl()[[1]]

form<-c(formA)                                                        

#User-defined model 
#form<-function(t,p,X){
#V<-p[1]
#Cl<-p[2]

#y<-X/V*(exp(-Cl/V*t))
#return(y)
#}

                                                                                                    