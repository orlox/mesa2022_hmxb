# This file was generated, do not modify it. # hide
#hideall
using ForwardDiff

# unitless potential, mass in units of M1, distance in units of a
function roche_potential(pos::Vector,q)
    x=pos[1]
    y=pos[2]
    z=pos[3]

    r = sqrt(x^2+y^2+z^2)
    rp = sqrt((1-x)^2+y^2+z^2)

    phi = -1/r - q*(1/rp-x) -(1+q)/2*(x^2+y^2)

    return phi
end

function find_Lpoint(q, point::Int)
    potential_func = x->roche_potential(x,q)

    #get gradient of pontential
    g = x -> ForwardDiff.gradient(potential_func, x)

    #solve Lagrangian points by bisection
    #set bounds for bisection solver
    xlow = 0
    xhigh = 0
    if point==1 #solve for L1
        xlow = 0
        xhigh = 1
    else
        #q<1 means M1>M2, x coordinate is centered at M1
        if q<1 #L2 is at positive X, L3 at negative
            if point==2
                xlow = 1
                xhigh = 3
            elseif point==3
                xlow = -2
                xhigh = 0
            else
                error("point must be 1,2 or 3")
            end 
        else #L2 is at negative X, L3 at positive
            if point==2
                xlow = -2
                xhigh = 0
            elseif point==3
                xlow = 1
                xhigh = 3
            else
                error("point must be 1,2 or 3")
            end 
        end
    end
    #print("$xlow,$xhigh\n")
    x_L = 0.5*(xlow+xhigh)
    n=0

    while abs(xlow-xhigh)>1e-10
        #evaluate x component of gradient. If positive, need to update xlow,
        #otherwise need to update xhigh
        grad = g([x_L,0,0])
        grad_x = grad[1]
        if grad_x > 0.0
            xlow = x_L
        else
            xhigh = x_L
        end

        #bisect limits
        x_L = 0.5*(xlow+xhigh)

        #print("$n $x_L $xlow $xhigh $grad_x \n")

        n+=1
        if n>10000
            error("Did not converge to solution")
        end

    end
    return x_L
end

# check things by plotting a few equipotentials in the X-Y plane
using Plots

q=0.5
x_L1 = find_Lpoint(q, 1)
x_L2 = find_Lpoint(q, 2)
x_L3 = find_Lpoint(q, 3)

potential_1D = (x)->roche_potential([x,0,0],q)
potential_2D = (x,y)->roche_potential([x,y,0],q)

phi_L1 = potential_1D(x_L1)
phi_L2 = potential_1D(x_L2)
phi_L3 = potential_1D(x_L3)

xvals = LinRange(-1.5,2,100)
plot(xvals, potential_1D.(xvals), legend=false)
xlabel!("\$x/a\$")
ylabel!("\$\\Phi/(GM_1/a)\$")
ylims!(-4,-2)
savefig(joinpath(@OUTPUT, "roche_1d.svg"))