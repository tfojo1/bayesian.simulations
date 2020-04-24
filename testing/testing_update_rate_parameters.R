

alpha=0.5
C=1
prior.iter=10

print(paste0('prior.iter=',prior.iter,', alpha=',alpha, ', C=',C))
all.next.iter = c(1,10,20,50,100,1000)

for (next.iter in all.next.iter)
{
#    rel.weight = round(sum(C/(prior.iter+1:next.iter)^alpha),2)
    gammas = C/(prior.iter+1:next.iter)^alpha
    p0 = prod(1-gammas)
    rel.weight = (1-p0)/p0

#    print(paste0(next.iter, '-->', rel.weight))
   # print(paste0(next.iter, '-->', round(rel.weight,5)))


    print(paste0(next.iter, '-->', round(1-p0,5)))

}
#    print(paste0('For first ', next.iter,
#                 ' iterations, weight is ',
#                 round(sum(C/(prior.iter+1:next.iter)^alpha),2),
#                 ' times the weight given to the initial value'))
#
