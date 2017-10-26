import PyPlot
const plt = PyPlot
plt.ion()
fig = plt.figure(1,figsize=(10,4))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
plt.tight_layout()

function plotR(ax,R,dims=[];reset=true)
	reset ? (ax[:cla](); ax[:plot]([0,1.2*R[:max]],[0,1.2*R[:max]],c="black",alpha=0.5)) : 0
	dims = typeof(dims) <: Int ? [dims] : (length(dims) > 0 ?  dims : R[:dims])
	# avg_birth = sum(R[:barcode][dim][:,1])/length(R[:barcode][dim][:,1])
	# avg_death = sum(R[:barcode][dim][:,2])/length(R[:barcode][dim][:,2])
    [ax[:scatter](R[:barcode][dim][:,1],R[:barcode][dim][:,2],s=5,alpha=0.5) for dim in dims]

end
