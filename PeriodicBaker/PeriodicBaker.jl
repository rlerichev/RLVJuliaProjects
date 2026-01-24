### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 4b8a3a21-889f-4189-a6c5-e06014e1e8e0
begin
	using Pkg
	urlsdd="https://github.com/Colectivo-SDD/SDD"
	Pkg.add(url=urlsdd*"Core.jl")
	Pkg.add(url=urlsdd*"Graphics.jl")
	Pkg.add(url=urlsdd*"Geometry.jl")
	Pkg.add(url=urlsdd*".jl", rev="rlv_makie")	
	#Pkg.add(url=urlsdd*"IFS.jl", rev="rlv_makie")
end

# ╔═╡ 7af6c311-4ff3-4f71-84f0-8d9ef2d591c8
using Colors, ColorSchemes, Images, ImageIO, FileIO

# ╔═╡ df10e677-8b09-473e-8b15-3ee063496f34
using Makie, GLMakie

# ╔═╡ 4294274f-8874-4840-a241-0cce4969d5dd
using SDDCore, SDDGraphics, SDDGeometry, SDD

# ╔═╡ d52d0e77-9a16-4163-b2ef-10779cd1cdff
include("../common/src/nbutils.jl")

# ╔═╡ afe66166-e9ca-11f0-bad5-59bd6959ad58
md"""
# Dynamics of $f_a$

### Function and iterates

$f_a(z)=\frac{-ze^z-1}{z-a}=\frac{ze^z+1}{a-z}$

Note that

$\lim_{t\rightarrow\infty}f_a(t)=\infty$

$\lim_{t\rightarrow-\infty}f_a(t)=0$

Second iterate

$f^2_a(z)=\frac{(\frac{ze^z+1}{a-z})e^\frac{ze^z+1}{a-z}+1}{a-\frac{ze^z+1}{a-z}}$

$=\frac{(ze^z+1)e^\frac{ze^z+1}{a-z}+a-z}{a(a-z)-ze^z-1}$
"""

# ╔═╡ f6b33556-cda7-4e17-8b1d-22801fb00c85
md"""
### Derivative and critical points

Derivative:

$f'_a(z)=\frac{(e^z+ze^z)(a-z)-(ze^z+1)(-1)}{(a-z)^2}$

$=\frac{ae^z+aze^z -ze^z-z^2e^z +ze^z+1}{(a-z)^2}$

$=\frac{(a+az-z^2)e^z+1}{(a-z)^2}$

Critical points:

$f'_a(z)=0 \implies (a+az-z^2)e^z+1=0$

$\implies -z^2+az+a+\frac{1}{e^z}=0$

Auxiliar function:

$g_a(z)=(a+az-z^2)e^z+1$

$g_a(z)=0 \implies f'_a(z)=0$

$g'_a(z)=(a-2z)e^z+(a+az-z^2)e^z=(2a+(a-2)z-z^2)e^z$
"""

# ╔═╡ b4c616d2-e720-465d-927d-e8534dc6cbec
function f(a::Number, z::Number)
	if isinf(z)
		if real(z)<0
			return 0
		else
			return Inf
		end
	elseif abs2(z-a)<0.00000001
		return Inf
	end
	(z*exp(z)+1)/(a-z)
end

# ╔═╡ 8e255a87-fe1c-4ab0-9330-1128cc317744
f´(a::Number, z::Number) = ((a+a*z-z^2)*exp(z)+1)/((a-z)^2)

# ╔═╡ e9a1d66a-abe1-4451-96ce-eccc6e887731
g(a::Number, z::Number) = (a+a*z-z^2)*exp(z)+1

# ╔═╡ 1464c2e2-aebf-488b-b1a5-15895edd0d05
g´(a::Number, z::Number) = (2a+(a-2)*z-z^2)*exp(z)

# ╔═╡ c7fa0c18-f71a-4afa-8541-5f00dbfad2b5
function createfa(a::Number)
	function f(z::Number)
		 (z*exp(z)+1)/(a-z)
	end
end

# ╔═╡ 7d58b085-2f4d-4f2f-a392-5074b98d7663
function createf´a(a::Number)
	function f(z::Number)
		 ((a+a*z-z^2)*exp(z)+1)/((a-z)^2)
	end
end

# ╔═╡ 8b152eed-60e1-4e0b-abfc-e38d2f06c1ef
md"""
### Preimages maps of $f'_a$ and $g_a$

The curves "ends" are the critical points, or the pole $a$.
"""

# ╔═╡ fecc94eb-e1f5-470a-8d5f-72d8c599f1f6
md"""
## Critical Points Approximation
"""

# ╔═╡ d6f5050a-0c02-481b-a376-bff33aeab725
function newtonraphsonC(f::Function, f´::Function, z0::Number=0.0; ε::Real=0.0000001, maxiterations=100)
	z = z0
	for n in 1:maxiterations
		if abs2(z) < ε
			return z
		end
		f´z = f´(z)
		if abs2(f´z) < ε
			return Inf
		end
		z = z - f(z)/f´z
	end
	z
end

# ╔═╡ e8e876c1-a196-406a-b444-a8618a208690
function compassC(f::Function, zC::Number, size::Real=1.0; ε::Real=0.0000001, maxiterations=16)

	function compassCRec(zC::Number, size::Real=1.0, level::Integer=16)
		
		mC = abs2(f(zC))
		
		if mC < ε || level == 0
			return zC
		end

		zS = zC - size*im
		zE = zC + size
		zN = zC + size*im
		zW = zC - size
		mS = abs2(f(zS))
		mE = abs2(f(zE))
		mN = abs2(f(zN))
		mW = abs2(f(zW))

		if mC <= mS && mC <= mE && mC <= mN && mC <= mW
			return compassCRec(zC, size/2, level-1)
		else
			if     mS <= mE && mS <= mN && mS <= mW
				return compassCRec(zS, size, level-1)
			elseif mE <= mN && mE <= mW && mE <= mS
				return compassCRec(zE, size, level-1)
			elseif mN <= mW && mN <= mS && mN <= mE
				return compassCRec(zN, size, level-1)
			elseif mW <= mS && mW <= mE && mW <= mN
				return compassCRec(zW, size, level-1)
			end
		end	

		return zC
	end

	compassCRec(zC, size, maxiterations)
end

# ╔═╡ 639e9512-7281-4bc7-8bfc-76f4409306c6
function compassC1(f::Function, zC::Number, size::Real=1.0; ε::Real=0.0000001, maxiterations=16)

	function compassCRec(zC::Number, size::Real=1.0, level::Integer=16)
		
		mC = abs(f(zC))
		
		if mC < ε || level == 0
			return zC
		end

		zS = zC - size*im
		zE = zC + size
		zN = zC + size*im
		zW = zC - size
		mS = abs(f(zS))
		mE = abs(f(zE))
		mN = abs(f(zN))
		mW = abs(f(zW))

		if mC <= mS && mC <= mE && mC <= mN && mC <= mW
			return compassCRec(zC, size/2, level-1)
		else
			if     mS <= mE && mS <= mN && mS <= mW
				return compassCRec(zS, size, level-1)
			elseif mE <= mN && mE <= mW && mE <= mS
				return compassCRec(zE, size, level-1)
			elseif mN <= mW && mN <= mS && mN <= mE
				return compassCRec(zN, size, level-1)
			elseif mW <= mS && mW <= mE && mW <= mN
				return compassCRec(zW, size, level-1)
			end
		end	

		return zC
	end

	compassCRec(zC, size, maxiterations)
end

# ╔═╡ b02bddf3-40c1-4cb9-bdee-7de17456bba1
g(-im,-im)

# ╔═╡ 208f363e-dc10-4307-aa27-3299ba2057a5
md"""
## Dynamic plane

Drawing

- Orbits:

$o(z,f)=\{z,f(z),f^2(z),\dots\}.$


- The trapped points set:

$\mathcal{K}(f)=\{z\in\mathbb{C}|\,\,o(z,f)\,\,\mathrm{is\, bounded}\}.$
"""

# ╔═╡ e0ff84af-11fd-4c81-8792-8c45e20984b2
function stops(f::Function, z0::Number, N::Integer=1; ε::Real=0.000001)
	z = z0
	fnz = z0
	for k in 1:N
		fnz = f(fnz)
	end
	if abs2(z-fnz) < ε
		return true
	end	
	false
end

# ╔═╡ ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
md"""
## Parameter plane
"""

# ╔═╡ 597ce818-bbd4-44b3-b985-80127d9feae7
let
	L(a,x)=a*x*(1-x)
	
	Npix = 600
	xmin,xmax,ymin,ymax = -0.1,4.1,-0.1,1.1
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	its1 = 10
	its2 = 10
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	for a in xs
		c0 = 0.01
		c1 = 0.5

		for n in 1:its1
			c0 = L(a,c0)
			c1 = L(a,c1)
		end

		points0 = [Point2f(a,real(c0))]		
		points1 = [Point2f(a,real(c1))]		
		for n in 1:its2
			c0 = L(a,c0)
			push!(points0, Point2f(a,real(c0)))
			c1 = L(a,c1)
			push!(points1, Point2f(a,real(c1)))
		end

		scatter!(ax, deepcopy(points0), markersize=2, color=:black)
		scatter!(ax, deepcopy(points1), markersize=2, color=:red)
		#scatter!(ax, points0, markersize=2, color=:black)
		#scatter!(ax, points1, markersize=2, color=:red)

		empty!(points0)
		empty!(points1)
	end

	fig
end

# ╔═╡ 7a0d58fd-7238-4f0c-accf-9e429ca4d325
md"""
## Software Julia packages
"""

# ╔═╡ a07e4983-d801-47ba-9531-16bdfb1f5f26
const Gr = SDDGraphics

# ╔═╡ c14eacb3-f173-4609-b750-413027a296ca
prismx = push!(RGBA.(deepcopy(colorschemes[:prism].colors)), RGBA(1,1,1,0.25))

# ╔═╡ fd700d57-b159-4a0e-87e6-6b7361254b60
vermeerx = pushfirst!(push!(RGBA.(Gr.reverse(:vermeer)), RGBA(0,0,0,0.9)),RGBA(1,1,1,0.1))

# ╔═╡ 57f7c4ba-d10a-4f0d-ae0d-b12f49a1c313
cuherx = Gr.reverse(:cubehelix) #pushfirst!(push!(RGBA.(Gr.reverse(:cubehelix)), RGBA(0,0,0,0.9)),RGBA(1,1,1,0.1))

# ╔═╡ 7087d5f4-86a4-4ed6-9970-3361ab584c79
let
	Npix = 720
	xmin,xmax,ymin,ymax = -9,6,-4,4
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = max(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	#maxits = 75
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,3Npix/4))
	ax3 = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			title=L"z\mapsto z",
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	#cf = ClassicDomainCF(cartesianshade=RGBA(0,0,0,0), cartesianexp=4)
	cf = AngleCF(0)
	preimages!(ax3, z->z, xs, ys,
			  coloringfunction = cf, colormap=cm,
			  interpolate=true, fxaa=true, ssao=true, depth_shift=1)

	fig
end

# ╔═╡ 559b6103-c4f9-4edb-9f30-fc51ca2184cf
let
	a = -1+4im

	Npix = 720
	xmin,xmax,ymin,ymax = -9,6,-4,24
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = max(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	#maxits = 75
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			title=L"z\mapsto((a+az-z^2)e^z+1)/(a-z)^2",
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)
	ax2 = Makie.Axis(fig[1,2], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			title=L"z\mapsto(a+az-z^2)e^z+1",
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	#cf = ClassicDomainCF(cartesianshade=RGBA(0,0,0,0), cartesianexp=4)
	cf = AngleCF(0)
	preimages!(ax, z->((a+a*z-z^2)*exp(z)+1)/((a-z)^2), xs, ys, iterations=1,
			  coloringfunction = cf, colormap=cm,
			  interpolate=true, fxaa=true, ssao=true, depth_shift=1)
	preimages!(ax2, z->((a+a*z-z^2)*exp(z)+1), xs, ys,
			  coloringfunction = cf, colormap=cm,
			  interpolate=true, fxaa=true, ssao=true, depth_shift=1)

	fig
end

# ╔═╡ f58338d6-a72a-4e23-8814-d3a88a5bb73f
let
	a = -0.5+0.25im

	Npix = 720
	xmin,xmax,ymin,ymax = -8,8,-6,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = max(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	#maxits = 75
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,3Npix/4))
	#=ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			title=L"z\mapsto((a+az-z^2)e^z+1)/(a-z)^2",
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)=#
	ax2 = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			title=L"z\mapsto(a+az-z^2)e^z+1",
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	#cf = ClassicDomainCF(cartesianshade=RGBA(0,0,0,0), cartesianexp=4)
	cf = AngleCF(0)
	#=preimages!(ax, z->((a+a*z-z^2)*exp(z)+1)/((a-z)^2), xs, ys, iterations=1,
			  coloringfunction = cf, colormap=cm,
			  interpolate=true, fxaa=true, ssao=true, depth_shift=1)=#
	preimages!(ax2, z->g(a,z), xs, ys,
			  coloringfunction = cf, colormap=cm,
			  interpolate=true, fxaa=true, ssao=true, depth_shift=1)

	ga = g(a,a)
	m = abs(ga)
	z0 = a+ga/m
	#c1 = compassC(z->g(a,z), real(a) >= 0 ? a+1 : g(a,a)/m, m/2) # a+1.5a/abs(a)
	c1 = compassC1(z->g(a,z), z0) #real(a) >= 0 ? a+2 : m, real(a) >= 0 ? 1 : 0.9m) # a+1.5a/abs(a)
	c1nr = newtonraphsonC(z->g(a,z), z->g´(a,z), z0)
	
	scatter!(ax2, [Point2f(real(a), imag(a))], color=:yellow, markersize=8)
	scatter!(ax2, [Point2f(real(z0), imag(z0))], color=:green, markersize=8)
	scatter!(ax2, [Point2f(real(c1), imag(c1))], color=:red, markersize=12)
	scatter!(ax2, [Point2f(real(c1nr), imag(c1nr))], color=:magenta, markersize=8)

	print(c1, ", ", c1nr)
	
	fig
end

# ╔═╡ 8bcefc7f-d7ba-4566-a3be-7400f2797964
let
	a = 1

	Npix = 600
	xmin,xmax,ymin,ymax = -6,6,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 60
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(3Npix/4,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	f = createfa(a)
	
	trappedpoints!(ax, f, xs, ys, maxiterations=maxits,
		hasescaped = z -> real(z) < -24, 
		colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#c1 = findcritic1b(l, m, maxiterations=20)
	#c1 = findcriticBif(0.001,2.999, l, m, maxiterations=20)
	#c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)
	#c3 = findcriticN(3, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	#=text!(ax, Point2f(real(c1),imag(c1)), text=L"c_1", color=:red, align=(:left,:center),
		  fontsize=0.5, markerspace=:data, offset=(0.25,0))
	text!(ax, Point2f(real(c2),imag(c2)), text=L"c_2", color=:blue, align=(:right,:center),
		  fontsize=0.5, markerspace=:data, offset=(-0.25,0))
	text!(ax, Point2f(real(c3),imag(c3)), text=L"c_3", color=:green, align=(:right,:center),
		  fontsize=0.5, markerspace=:data, offset=(-0.25,0))
	=#

	#c1 = compassC(z->g(a,z), real(a) >= 0 ? a+1 : g(a,a)/abs(g(a,a)), abs(g(a,a))/2)	
	ga = g(a,a)
	z0 = a + ga/abs(ga)
	c1 = newtonraphsonC(z->g(a,z), z->g´(a,z), z0)
	
	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.75)], linewidth=0.5, markersize=4, iterations=6)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	#scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	#scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	#scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	#text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
	#	  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m025.jpg", fig)
	
	fig
end

# ╔═╡ 0947a4b6-0309-4fea-bcab-7a3f0d036974
cuhex = :cubehelix #pushfirst!(push!(RGBA.(colorschemes[:cubehelix].colors), RGBA(1,1,1,0.1)),RGBA(0,0,0,0.9))

# ╔═╡ 7e516588-5ad2-498b-89a6-91a376600ab2
cuherx_dark = pushfirst!(Gr.reverse(:cubehelix),RGB(0,0,0))

# ╔═╡ 8a8619e8-0236-4c37-bd09-65d82110fc2a
md"""
## Critical points: Intersection of curves

(Back up reference from the case $f_{\lambda,\mu}(z)=\lambda e^z + \mu/z$)
"""

# ╔═╡ dcc54357-7342-4978-a0e3-ae25b3161158
f(l::Number, m::Number, z::Number) = l*exp(z)+m/z

# ╔═╡ 97f0ae6e-02b2-4f5f-a0d6-83c1bbac8666
f(1,-Inf)

# ╔═╡ ab535ab3-bb51-496d-b984-aef9a1b991d1
let
	Npix = 400
	xmin,xmax,ymin,ymax = -2,4.5,-1.25,1.25
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 100
	cm = cuherx #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,9Npix/16))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.25:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(a,z) -> f(a, z), xs, ys,
		#seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		seed = -0.1,
		#hasescaped = (a,z) -> real(z)>36,
		#hasescaped = (a,z) -> abs2(z-f(a,z))<0.001,
		hasescaped = (a,z) -> real(z)>64 || abs2(z-f(a,z))<0.0001,
		#hasescaped = (a,z) -> abs(real(z))>36 || abs2(z-f(a,f(a,z)))<0.0001,
		#hasescaped = (a,z) -> abs2(z-f(a,f(a,z)))<0.0001,
		#hasescaped = (a,z) -> abs2(z-f(a,f(a,f(a,z))))<0.0001,
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)
	#scatter!(ax, [Point2f(-20,0.25)], markersize=12, color=:red)
	#text!(Point2f(-20,0.25), text="(-20,0.25)", color=:red, align = (:left,:bottom), fontsize=0.5,
#		markerspace= :data)	

	Colorbar(fig[2, 1], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = false, flipaxis = false)

	#save("Sienra_Leaf_0.jpg", fig)
	
	fig
end

# ╔═╡ d5872959-a8a4-443a-909a-265eddb9b5a1
let
	Npix = 1200
	xmin,xmax,ymin,ymax = -2,4.5,-1.25,1.25
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,9Npix/16))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.25:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(a,z) -> f(a, z), xs, ys,
		seed = a -> newtonraphsonC(z->g(a,z), z->g´(a,z), a+g(a,a)/abs(g(a,a))), 
		#seed = a -> compassC(z->g(a,z), real(a) >= 0 ? a+1 : g(a,a)/abs(g(a,a)), abs(g(a,a))/2),
		#seed = a -> compassC(z->g(a,z), real(a) >= 0 ? a+2 : abs(g(a,a)), real(a) >= 0 ? 1 : 0.9abs(g(a,a))),
		#hasescaped = (a,z) -> real(z)>36,
		#hasescaped = (a,z) -> abs2(z-f(a,z))<0.001,
		hasescaped = (a,z) -> real(z)>144 || abs2(z-f(a,z))<0.0001,
		#hasescaped = (a,z) -> abs(real(z))>36 || abs2(z-f(a,f(a,z)))<0.0001,
		#hasescaped = (a,z) -> abs2(z-f(a,f(a,z)))<0.0001,
		#hasescaped = (a,z) -> abs2(z-f(a,f(a,f(a,z))))<0.0001,
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)
	#scatter!(ax, [Point2f(-20,0.25)], markersize=12, color=:red)
	#text!(Point2f(-20,0.25), text="(-20,0.25)", color=:red, align = (:left,:bottom), fontsize=0.5,
#		markerspace= :data)	

	Colorbar(fig[2, 1], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = false, flipaxis = false)

	save("SienraLeaf_critic.jpg", fig)
	
	fig
end

# ╔═╡ 683e2ae4-07cd-48b2-839f-717896f3e8a1
let
	Npix = 800
	xmin,xmax,ymin,ymax = -3,5,-4,4
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = 0.002 #min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	its1 = 20
	its2 = 20
	cm = cuherx #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.5:ymax, ygridcolor=gridc)

	for a in xs
		c0 = 0
		ga = g(a,a)
		z0 = a + ga/abs(ga)
		c1 = newtonraphsonC(z->g(a,z), z->g´(a,z), z0)

		for n in 1:its1
			c1 = f(a,c1)
			c0 = f(a,c0)
		end

		points0 = [Point2f(a,real(c0))]		
		points1 = [Point2f(a,real(c1))]		
		for n in 1:its2
			c0 = f(a,c0)
			push!(points0, Point2f(a,real(c0)))
			c1 = f(a,c1)
			push!(points1, Point2f(a,real(c1)))
		end

		#scatter!(ax, points0, markersize=4, color=:black)
		#scatter!(ax, points1, markersize=4, color=:red)
		scatter!(ax, deepcopy(points0), markersize=2, color=:black)
		scatter!(ax, deepcopy(points1), markersize=2, color=:red)
		
		empty!(points0)
		empty!(points1)
	end

	fig
end

# ╔═╡ a87d4962-eb05-45e4-9cc5-5367dd13b8d2
let
	Npix = 800
	xmin,xmax,ymin,ymax = -1.75,0,-2,1
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = 0.001 #min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	its1 = 40
	its2 = 20
	cm = cuherx #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.5:ymax, ygridcolor=gridc)

	for a in xs
		c0 = 0
		ga = g(a,a)
		z0 = a + ga/abs(ga)
		c1 = newtonraphsonC(z->g(a,z), z->g´(a,z), z0)

		for n in 1:its1
			c1 = f(a,c1)
			c0 = f(a,c0)
		end

		points0 = [Point2f(a,real(c0))]		
		points1 = [Point2f(a,real(c1))]		
		for n in 1:its2
			c0 = f(a,c0)
			push!(points0, Point2f(a,real(c0)))
			c1 = f(a,c1)
			push!(points1, Point2f(a,real(c1)))
		end

		#scatter!(ax, points0, markersize=4, color=:black)
		#scatter!(ax, points1, markersize=4, color=:red)
		scatter!(ax, deepcopy(points0), markersize=2, color=:black)
		scatter!(ax, deepcopy(points1), markersize=2, color=:red)
		
		empty!(points0)
		empty!(points1)
	end

	fig
end

# ╔═╡ 9acb88d5-f6a1-403a-ab7f-41501af82487
let
	Npix = 800
	xmin,xmax,ymin,ymax = 3.5,4.5,-0.5,1.5
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = 0.001 #min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	its1 = 40
	its2 = 20
	cm = cuherx #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.5:ymax, ygridcolor=gridc)

	for a in xs
		c0 = 0
		ga = g(a,a)
		z0 = a + ga/abs(ga)
		c1 = newtonraphsonC(z->g(a,z), z->g´(a,z), z0)

		for n in 1:its1
			c1 = f(a,c1)
			c0 = f(a,c0)
		end

		points0 = [Point2f(a,real(c0))]		
		points1 = [Point2f(a,real(c1))]		
		for n in 1:its2
			c0 = f(a,c0)
			push!(points0, Point2f(a,real(c0)))
			c1 = f(a,c1)
			push!(points1, Point2f(a,real(c1)))
		end

		#scatter!(ax, points0, markersize=4, color=:black)
		#scatter!(ax, points1, markersize=4, color=:red)
		scatter!(ax, deepcopy(points0), markersize=2, color=:black)
		scatter!(ax, deepcopy(points1), markersize=2, color=:red)
		
		empty!(points0)
		empty!(points1)
	end

	fig
end

# ╔═╡ 04d2ed2f-de6c-4e5a-a5de-c848b6de8a87
f2(l::Number, m::Number, z::Number) = l*exp(l*exp(z)+m/z)+m*z/(l*z*exp(z)+m)

# ╔═╡ 729aaaf9-28f3-4274-9a26-12fc6f80f074
f´(l::Number, m::Number, z::Number) = l*exp(z)-m/(z^2)

# ╔═╡ 4662f2ab-8b41-4211-98e0-58778d61d007
f´´(l::Number, m::Number, z::Number) = l*exp(z)+2m/(z^3)

# ╔═╡ f09e37b5-b0af-46c3-bf25-51017bf48ffa
function createflm(l::Number, m::Number)
	function f(z::Number)
		 l*exp(z)+m/z
	end
end

# ╔═╡ a2cf3c8b-c5c8-4df5-951e-62a7954e1083
function createflm2(l::Number, m::Number)
	function f(z::Number)
		lez = l*exp(z)
		l*exp(lez+m/z)+m*z/(lez*z+m)
	end
end

# ╔═╡ 9bceda1a-3a79-4e35-95b2-2572d5577993
function createf´lm(l::Number, m::Number)
	function f(z::Number)
		 l*exp(z)-m/(z^2)
	end
end

# ╔═╡ cada7f65-ccb1-4ec3-98e6-2370aac392cc
function createf´´lm(l::Number, m::Number)
	function f(z::Number)
		 l*exp(z)+2m/(z^3)
	end
end

# ╔═╡ 502034db-8faa-41ae-92be-c150305a1105
function cy(x::Real, l::Real, m::Real)
	Δ = abs(m/l)/exp(x)-x^2
	if Δ < 0.0
		return 0.0
	end
	sqrt(Δ)
end

# ╔═╡ 7e62773c-550b-483e-bf86-efed99abbbc6
cx(y::Real, l::Real, m::Real) = -y/(tan((y+pi)/2))

# ╔═╡ 3031e81a-5cd4-47f5-8d24-0d407ea29327
cxC(y::Real, l::Number, m::Number) = y/(tan((angle(m/l)-y)/2))

# ╔═╡ 2f247a63-5943-49e7-a71e-f6d94af13181
md"""
### Finding critical points...
"""

# ╔═╡ 75927946-3b6a-4a06-ba52-669066a434f5
function findcritic1(x0::Real, l::Real, m::Real, 
	ε::Real=0.0000001, maxiterations=100)
	
	#x,y = cx(y0,l,m),y0
	x,y = x0, cy(x0,l,m)
	crit = complex(x,y)
	n = 0
	while abs2(f´(l,m,crit)) > ε && n < maxiterations
		x = cx(y,l,m)
		y = cy(x,l,m)
		crit = complex(x,y)
		n = n+1
	end

	#print(n)
	return crit
end

# ╔═╡ 6dce489a-ec3d-4907-b2b6-01e9043fe394
function findcritic1b(l::Real, m::Real; ε::Real=0.0000001, maxiterations=100)
	yl = 0.001
	yr = 2.999
	#yms = []
	crit = complex(cx(yl,l,m), yl)
	for n in 1:maxiterations
		ym = (yl+yr)/2
		#push!(yms,[ym,n])
		xm = cx(ym,l,m)
		Δ = abs(m/l)*exp(xm)-xm^2
		if Δ < 0
			yr = ym
		else
			x = cx(ym,l,m)
			y = cy(x,l,m)
			crit = complex(x,y)
			if abs2(f´(l,m,crit)) < ε
				return crit
			end
			if y > ym
				yl = ym
			else
				yr = ym
			end
		end
	end
	#print(yms)
	return crit
end

# ╔═╡ 1e61565e-072a-4d53-a573-41f91f191074
function findcriticBif(yleft::Real, yright::Real, l::Real, m::Real; ε::Real=0.0000001, maxiterations=100)
	yl = yleft
	yr = yright
	ym = (yl+yr)/2
	crit = complex(cx(ym,l,m), ym)
	for n in 1:maxiterations
		ym = (yl+yr)/2
		xm = cx(ym,l,m)
		Δ = abs(m/l)*exp(xm)-xm^2
		if Δ < 0
			yr = ym
		else
			x = cx(ym,l,m)
			y = cy(x,l,m)
			#crit = complex(x,y)
			#if abs2(f´(l,m,crit)) < ε
			if abs2(y-ym) < ε
				#print(n)
				return complex(x,y)
			end
			if y > ym
				yl = ym
			else
				yr = ym
			end
		end
	end
	x = cx(ym,l,m)
	y = cy(x,l,m)
	return crit
end

# ╔═╡ 978cc14e-3383-450b-a46f-0376baa84b8c
function newtonraphsonR(f::Function, f´::Function, x0::Real=0.0; ε::Real=0.0000001, maxiterations=100)
	x = x0
	for n in 1:maxiterations
		if abs(x) < ε
			return x
		end
		f´x = f´(x)
		if abs(f´x) < ε
			return Inf
		end
		x = x - f(x)/f´x
	end
	x
end

# ╔═╡ Cell order:
# ╟─afe66166-e9ca-11f0-bad5-59bd6959ad58
# ╟─f6b33556-cda7-4e17-8b1d-22801fb00c85
# ╠═b4c616d2-e720-465d-927d-e8534dc6cbec
# ╠═97f0ae6e-02b2-4f5f-a0d6-83c1bbac8666
# ╠═8e255a87-fe1c-4ab0-9330-1128cc317744
# ╠═e9a1d66a-abe1-4451-96ce-eccc6e887731
# ╠═1464c2e2-aebf-488b-b1a5-15895edd0d05
# ╠═c7fa0c18-f71a-4afa-8541-5f00dbfad2b5
# ╠═7d58b085-2f4d-4f2f-a392-5074b98d7663
# ╟─8b152eed-60e1-4e0b-abfc-e38d2f06c1ef
# ╟─7087d5f4-86a4-4ed6-9970-3361ab584c79
# ╠═559b6103-c4f9-4edb-9f30-fc51ca2184cf
# ╟─fecc94eb-e1f5-470a-8d5f-72d8c599f1f6
# ╠═d6f5050a-0c02-481b-a376-bff33aeab725
# ╟─e8e876c1-a196-406a-b444-a8618a208690
# ╟─639e9512-7281-4bc7-8bfc-76f4409306c6
# ╠═f58338d6-a72a-4e23-8814-d3a88a5bb73f
# ╠═b02bddf3-40c1-4cb9-bdee-7de17456bba1
# ╟─208f363e-dc10-4307-aa27-3299ba2057a5
# ╠═8bcefc7f-d7ba-4566-a3be-7400f2797964
# ╠═e0ff84af-11fd-4c81-8792-8c45e20984b2
# ╟─ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
# ╠═ab535ab3-bb51-496d-b984-aef9a1b991d1
# ╠═d5872959-a8a4-443a-909a-265eddb9b5a1
# ╠═683e2ae4-07cd-48b2-839f-717896f3e8a1
# ╠═a87d4962-eb05-45e4-9cc5-5367dd13b8d2
# ╠═9acb88d5-f6a1-403a-ab7f-41501af82487
# ╠═597ce818-bbd4-44b3-b985-80127d9feae7
# ╠═d52d0e77-9a16-4163-b2ef-10779cd1cdff
# ╠═7a0d58fd-7238-4f0c-accf-9e429ca4d325
# ╠═7af6c311-4ff3-4f71-84f0-8d9ef2d591c8
# ╠═df10e677-8b09-473e-8b15-3ee063496f34
# ╠═4b8a3a21-889f-4189-a6c5-e06014e1e8e0
# ╠═4294274f-8874-4840-a241-0cce4969d5dd
# ╠═a07e4983-d801-47ba-9531-16bdfb1f5f26
# ╠═c14eacb3-f173-4609-b750-413027a296ca
# ╠═fd700d57-b159-4a0e-87e6-6b7361254b60
# ╠═57f7c4ba-d10a-4f0d-ae0d-b12f49a1c313
# ╠═0947a4b6-0309-4fea-bcab-7a3f0d036974
# ╠═7e516588-5ad2-498b-89a6-91a376600ab2
# ╠═8a8619e8-0236-4c37-bd09-65d82110fc2a
# ╠═dcc54357-7342-4978-a0e3-ae25b3161158
# ╠═04d2ed2f-de6c-4e5a-a5de-c848b6de8a87
# ╠═729aaaf9-28f3-4274-9a26-12fc6f80f074
# ╠═4662f2ab-8b41-4211-98e0-58778d61d007
# ╠═f09e37b5-b0af-46c3-bf25-51017bf48ffa
# ╠═a2cf3c8b-c5c8-4df5-951e-62a7954e1083
# ╠═9bceda1a-3a79-4e35-95b2-2572d5577993
# ╠═cada7f65-ccb1-4ec3-98e6-2370aac392cc
# ╠═502034db-8faa-41ae-92be-c150305a1105
# ╠═7e62773c-550b-483e-bf86-efed99abbbc6
# ╠═3031e81a-5cd4-47f5-8d24-0d407ea29327
# ╟─2f247a63-5943-49e7-a71e-f6d94af13181
# ╠═75927946-3b6a-4a06-ba52-669066a434f5
# ╠═6dce489a-ec3d-4907-b2b6-01e9043fe394
# ╠═1e61565e-072a-4d53-a573-41f91f191074
# ╠═978cc14e-3383-450b-a46f-0376baa84b8c
