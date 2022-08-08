# This file was generated, do not modify it. # hide
#hideall

# unitless potential, mass in units of M1, distance in units of a
function a_div_ai(qf,qi)
    return (qi/qf)^(2)*((1+qi)/(1+qf))*exp(2*(qf-qi))
end

# check things by plotting a few equipotentials in the X-Y plane
using Plots

Mdonor = 30
Mfinal = 15


xvals = LinRange(15,30,100)
plot()
for Ma in [4,5,6]
    qvals = xvals./Ma
    qi = 30/Ma
    plot!(xvals, map(qf->log10(a_div_ai(qf,qi)), qvals), label=string(Ma)*"\$M_\\odot\$", legend=:bottomright)
end
xlabel!("\$M_\\mathrm{donor}\\;[M_{\\odot}]\$")
ylabel!("\$\\log(a/a_\\mathrm{i})\$")
savefig(joinpath(@OUTPUT, "shrinkage.svg"))