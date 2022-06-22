@def title = "MESA summer school 2022: Modeling of HMXBs with MESA"
@def tags = ["syntax", "code"]

```julia:plotsetup
#hideall
using Plots

size_inches = (3.38, 2.535)
size_pt = 1.5*72 .* size_inches

theme(:default;
    linewidth=3,
    fontfamily = "Computer Modern",
    size=size_pt,

    colorbar_tickfontsize = 1,
    xtickfontsize=9,
    ytickfontsize=9,
    xlabelfontsize=10,
    ylabelfontsize=10,
    
    minorticks=4,
    grid = false,
    
    framestyle=:box
)
```

# Welcome

\tableofcontents <!-- you can use \toc as well -->

## One part

$$F(x)=\int_0^1 f(x)dx$$

Some text with inline math $a=1+1$

A fortran snippet
```fortran
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
```

An inlist snippet
```fortran
&binary_job

   inlist_names(1) = 'inlist1' 
   inlist_names(2) = 'inlist2'

   evolve_both_stars = .true.
   ! rlof is turned on once star definitely settles into ZAMS
   change_initial_ignore_rlof_flag = .true.
   new_ignore_rlof_flag = .true.

/ ! end of binary_job namelist

&binary_controls
   do_tidal_sync = .true.
   do_j_accretion = .true.
   ! be 100% sure MB is always off
   do_jdot_mb = .false.
   do_jdot_missing_wind = .true.

   terminate_if_initial_overflow = .false.

   mdot_scheme = "contact"
   report_rlo_solver_progress = .true.

   ! initial conditions specified in extra inlist
   ! the initial masses and period are given there
   read_extra_binary_controls_inlist1 = .true.
   extra_binary_controls_inlist1_name = "inlist_extra"

   terminal_interval = 10

   ! timestep controls, these are relaxed for the test suite, recommended to use the commented
   ! values or lower for serious work
   fr = 0.1
   !fr = 0.02
   fr_limit = 1d-2
   fm = 0.1
   fm_limit = 1d-1
   varcontrol_case_a = 1d-3
   !varcontrol_case_a = 4d-4
   varcontrol_case_b = 1d-3
   !varcontrol_case_b = 5d-4
   varcontrol_ms = 1d-3
   !varcontrol_ms = 4d-4
   varcontrol_post_ms = 1d-3
   !varcontrol_post_ms = 5d-4
   dt_softening_factor = 0.4

   limit_retention_by_mdot_edd = .false.
   implicit_scheme_tolerance = 2d-4
   max_tries_to_achieve = 400
   min_change_factor = 1.1d0
   max_change_factor = 1.5d0
   initial_change_factor = 1.2d0
   change_factor_fraction = 0.8d0
   implicit_lambda = 0.4d0
   roche_min_mdot = 1d-10

   sync_mode_1 = "Uniform"
   sync_type_1 = "Hut_rad"
   Ftid_1 = 1
   sync_mode_2 = "Uniform"
   sync_type_2 = "Hut_rad"
   Ftid_2 = 1
   do_initial_orbit_sync_1 = .true.
   do_initial_orbit_sync_2 = .true.

   min_mdot_for_implicit = 1d-10
   roche_min_mdot = 1d-12
   accretor_overflow_terminate = 100.0d0
         
/ ! end of binary_controls namelist
```

## Another part

```julia:plot_roche
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
phi_L4 = potential_2D(0.5,sqrt(1-0.5^2))

equipotentials = [2*phi_L1, phi_L1, 0.5*(phi_L1+phi_L2),phi_L2,0.5*(phi_L2+phi_L3),phi_L3,0.5*(phi_L3+phi_L4),0.2*phi_L3+0.8*phi_L4]

xvals = LinRange(-1,2,500)
yvals = LinRange(-1,1,500)
contour(xvals, yvals, potential_2D, levels=equipotentials, color=:black,colorbar=false, linewidth=1)
scatter!([x_L1,x_L2,x_L3,0.5,0.5],[0,0,0,sqrt(1-0.5^2),-sqrt(1-0.5^2)],markersize=8,label=nothing)
xlabel!("x/a")
ylabel!("y/a")
savefig(joinpath(@OUTPUT, "roche.svg"))

```
\figenv{the caption}{/assets/index/code/output/roche.svg}{width:100%;border: 1px solid gray;}

More text