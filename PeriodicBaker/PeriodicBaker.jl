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
f(a::Number, z::Number) = (z*exp(z)+1)/(a-z)

# ╔═╡ 8e255a87-fe1c-4ab0-9330-1128cc317744
f´(a::Number, z::Number) = ((a+a*z-z^2)*exp(z)+1)/((a-z)^2)

# ╔═╡ e9a1d66a-abe1-4451-96ce-eccc6e887731
g(a::Number, z::Number) = (a+a*z-z^2)*exp(z)+1

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

# ╔═╡ dcc54357-7342-4978-a0e3-ae25b3161158
f(l::Number, m::Number, z::Number) = l*exp(z)+m/z

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

# ╔═╡ 8a8619e8-0236-4c37-bd09-65d82110fc2a
md"""
## Critical points: Intersection of curves ???
"""

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

# ╔═╡ fdbde74c-6ede-422c-955f-de5e02ec8b22
let
	l,m = -10,2
	xmin,xmax,ymin,ymax = -1,1,0,2
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))
	#vlines!(ax, [m/l], color=:black, linestyle=:dash)
	
	xs = xmin:0.001:xmax
	lines!(ax, xs, cy.(xs,l,m), color=:cadetblue, linewidth=4)
	ys = ymin:0.001:ymax
	lines!(ax, cx.(ys,l,m), ys, color=:orange, linewidth=4)

	#c1 = findcritic1(-0.25m/l, l, m)
	#c1 = compassC(z -> f(l, m, z), -0.25m/l + cy(-0.25m/l, l, m)*im, 0.01)
	#c1 = findcritic1b(l, m, maxiterations=20, ε=0.0001)
	c1 = findcriticBif(0.001, 2.999, l, m) #, maxiterations=100, ε=0.00000001)
	#print(c1, ", ", abs2(f´(l,m,c1)))
	scatter!(ax, [Point2f(real(c1), imag(c1))], color=:red, markersize=16)
	
	fig
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

# ╔═╡ 781de362-0f5f-44e5-87cc-d77345617ca7
#compassC(z->(z-3)*(z-im)*(z+3+im), -3-0.9im, 1)

# ╔═╡ d2a17754-e606-48e9-82a1-a350091538d2
cΔy(x::Real, l::Real, m::Real) = abs(m/l)/exp(x)-x^2

# ╔═╡ 0bce24f6-9278-464e-9d0c-0af471089dec
cΔy´(x::Real, l::Real, m::Real) = -abs(m/l)/exp(x)-2x

# ╔═╡ e3130dbc-4b0d-4fa0-80a2-2980fd875d0b
newtonraphsonR(x->cΔy(x,-20,0.25), x->cΔy´(x,-20,0.25), -60)

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

# ╔═╡ 6b7870f5-a79d-416c-9a2d-6bd8a38e5fc5
function findcriticN(N::Integer, l::Real, m::Real; ε::Real=0.0000001, maxiterations=100, maxiterationscompass=16)
	x = newtonraphsonR(x->cΔy(x,l,m), x->cΔy´(x,l,m), 0.5*l/m, ε=ε, maxiterations=maxiterations)
	compassC(z->f´(l,m,z), complex(x,(2N-3)*pi + pi/2), 1, ε=ε, maxiterations=maxiterationscompass)
end

# ╔═╡ 4cc91117-7c7b-4e70-820c-aff85860ab02
let
	l,m = -2,0.25
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(-12,2,0,12))
	xs = -12:0.02:2
	lines!(ax, xs, cy.(xs,l,m), color=:cadetblue, linewidth=4)
	ys = 0:0.02:12
	lines!(ax, cx.(ys,l,m), ys, color=:orange, linewidth=4)
	hlines!(ax, [pi,2pi,3pi], color=:black, linestyle=:dash)

	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=24, ε=0.00000001)
	#c2 = findcriticBif(3.2, 5.8, l, m, maxiterations=400, ε=0.0000000001)
	scatter!(ax, [Point2f(real(c2), imag(c2))], color=:red, markersize=16)
	
	fig
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

# ╔═╡ e410e775-a417-40a2-8dbf-4c9c5050eef7
md"""
$(\lambda,\mu)=(-20,1/4)$
"""

# ╔═╡ 7273ed9e-5dfb-475b-b2a1-a65bfb64ab87
md"""
$(\lambda,\mu)=(-1,1/4)$

Attracting fixed.
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

# ╔═╡ cdf5ac1b-df7f-4270-8c67-1340e20861b4
md"""
In black, the "Mandelbrot set" of $f_{\lambda,\mu}$:

$\mathcal{M}_1=\{(\lambda,\mu)\in\mathbb{R}^2|\,\,o(c_1,f_{\lambda,\mu})\,\,\mathrm{is\,bounded}\}.$

In colors, using the escape time:

$\mathbb{R}^2-\mathcal{M}_1.$
"""

# ╔═╡ 0fd2c626-52e9-4f8c-8eac-64a2047d9729
md"""
Reference: The Mandelbrot set

$\mathcal{M}=\{c\in\mathbb{C}|\,\,o(0,q_c)\,\,\mathrm{is\,bounded}\}=\{z\in\mathbb{C}|\,\,q_c^n(0)\nrightarrow\infty\},$

where $q_c(z)=z^2+c$.
"""

# ╔═╡ d91916f8-37a6-4494-ac61-6df7aba6a0c1
#=let
	xmin,xmax,ymin,ymax = -2,0,0,2
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:0.25:xmax, yticks=ymin:0.25:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticN(2, real(lm), imag(lm), maxiterations=24, ε=0.000001),
		hasescaped=(c,z)->real(z)<-16, 
		maxiterations=80, colormap=Gr.reverse(:vangogh)
	)

	fig
end=#

# ╔═╡ 6e0bf940-1d94-4efa-a0a0-6558ba5d67fe
#=let
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:1:xmax, yticks=ymin:1:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.00001), # || real(z)<-16, 
		maxiterations=124, colormap=:vangogh
	)

	fig
end=#

# ╔═╡ 65365633-6555-4f3d-945b-66275dcf8c7b
md"""
Reference: Mandelbrot set.
"""

# ╔═╡ ffbacbc1-1b56-4b1f-93f1-0dce2f44bc9f
md"""
Let $z_0=z_0(\lambda,\mu)\in\mathbb{C}$ be the convergence point of the sequence $f_{\lambda,\mu}^{2n}(c_1)$, such that $|(f_{\lambda,\mu}^{2n})'(z_0)|<1$ or $z_0=0$.

In colors (coloring with $n\mapsto|f_{\lambda,\mu}^n(c_1)-z_0|<\varepsilon$):

$\mathcal{M}_2=\{(\lambda,\mu)\in\mathbb{R}^2|\,\, f_{\lambda,\mu}^{2n}(c_1)\rightarrow z_0\}.$

In black:

$\mathbb{R}^2-\mathcal{M}_2.$

In red, the line of parameters where $f_{\lambda,\mu}$ has fixed parabolic points with multiplier $-1$:

$\mathcal{L}_{-1}=\{(\lambda,\mu)\in\mathbb{R}^2|\,\,\exists z:\,\,f_{\lambda,\mu}(z)=z,\,\,f'_{\lambda,\mu}(z)=-1\}.$
"""

# ╔═╡ 249dda01-d980-4cc6-86f0-1df2a31dd7d6
md"""
Reference: in colors the Mandelbrot set's bulbs of period $1$ and $2$, union its exterior:

$\{z\in\mathbb{C}|\,\,q_c^{2n}(0)\rightarrow z_c\,\,\mathrm{or}\,\,q_c^{2n}(0)\rightarrow\infty\},$

where $q_c(z)=z^2+c$ and $z_c$ is an attracting fixed point of $q^2_c$.
"""

# ╔═╡ 07765e47-4d63-4df8-ab8f-906d71a75338
mandelbrot((c,z)->z^2+c, -2.01:0.002:0.51,-1.16:0.002:1.16, seed=0, colormap=:cubehelix, #vermeerx,
	hasescaped = (c,z) -> stops(z->z^2+c, z, 2, ε=0.00001) || abs2(z)>144, maxiterations=128,
	depth_shift=1, interpolate=true, fxaa=true, ssao=true, # transparency=false, alpha=0.9, 
	axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16), xticks=-2:0.25:0.5, yticks=-1:0.25:1, xgridcolor=RGBA(0.25,0.25,0.25,0.25), ygridcolor=RGBA(0.25,0.25,0.25,0.25), backgroundcolor=RGBA(1,1,1,0)))

# ╔═╡ 2b61308c-3efa-47a9-88d8-0ea870ca3d77
#=mandelbrot((c,z)->z^2+c, -2.01:0.002:0.51,-1.16:0.002:1.16, seed=0, colormap=cuherx, #vermeerx,
	hasescaped = (c,z) -> stops(z->z^2+c, z, 2, ε=0.00001) || abs2(z)>144, maxiterations=128, axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16)))
=#

# ╔═╡ 855b2724-3129-4166-b170-99f6ae6cbc37
md"""
Let $z_0=z_0(\lambda,\mu)\in\mathbb{C}$ be the convergence point of the sequence $f_{\lambda,\mu}^{2n}(c_1)$, such that $|(f_{\lambda,\mu}^{2n})'(z_0)|<1$ and $z_0\neq0$.

In colors (coloring with $n\mapsto|f_{\lambda,\mu}^n(c_1)-z_0|<\varepsilon$):

$\mathcal{M}_{2,z_0\neq0}=\{(\lambda,\mu)\in\mathbb{R}^2|\,\, f_{\lambda,\mu}^{2n}(c_1)\rightarrow z_0\}\subset\mathcal{M}_2.$

In black:

$\mathbb{R}^2-\mathcal{M}_{2,z_0\neq0}.$

In blue, the line $\mathcal{L}_{-1}$.
"""

# ╔═╡ ab526193-7dfd-4e82-a397-beb7e2fd6c93
md"""
Reference: Bulbs of period1 and 2 of the Mandelbrot set.
"""

# ╔═╡ e23f79d7-5899-409b-b721-4bb06aac37b4
mandelbrot((c,z)->z^2+c, -2.01:0.002:0.51,-1.16:0.002:1.16, seed=0, colormap=:cubehelix, #vermeerx,
	hasescaped = (c,z) -> stops(z->z^2+c, z, 2, ε=0.00001), maxiterations=128, #  || abs2(z)>144
	depth_shift=1, interpolate=true, fxaa=true, ssao=true, # transparency=false, alpha=0.9, 
	axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16), xticks=-2:0.25:0.5, yticks=-1:0.25:1, xgridcolor=RGBA(0.25,0.25,0.25,0.25), ygridcolor=RGBA(0.25,0.25,0.25,0.25), backgroundcolor=RGBA(1,1,1,0)))

# ╔═╡ 40f1cf47-2b5c-4d9e-bc8b-10ca781dffcc
md"""
Bifurcations, KLM paper.
"""

# ╔═╡ 05efe4d3-4126-449a-ac48-1799c4c2f28a
#=let
	xmin,xmax,ymin,ymax = -5,-3,4,5
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:0.25:xmax, yticks=ymin:0.25:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.00001), # || abs2(z)>100, 
		maxiterations=320, colormap=vermeerx #:vangogh
	)

	fig
end=#

# ╔═╡ de1c003f-d951-47cf-acc3-beaea021298b
# "Vortex"
#=let
	xmin,xmax,ymin,ymax = -2,-1,4,5
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:1:xmax, yticks=ymin:1:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.00001), # || real(z)<-16, 
		maxiterations=124, colormap=:vangogh
	)

	fig
end=#

# ╔═╡ 2c9f1a59-ad76-4476-8d32-7e158dcd2541
#=let
	xmin,xmax,ymin,ymax = -5,0,0,3
	Δx,Δy = (xmax-xmin)/800, (ymax-ymin)/800
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(800,600))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:0.5:xmax, yticks=ymin:0.5:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.000001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025), # || real(z)<-16, 
		maxiterations=200, colormap=vermeerx #Gr.reverse(:vermeer) #prismx #:vangogh
	)

	lines!(ax, [-ℯ,0], [0,1], color=:red, linewidth=1.5)

	fig
end=#

# ╔═╡ ee95a23b-177b-474d-a517-5cc9204831dd


# ╔═╡ fab68d9a-311d-451f-9289-caaec83569dc
marcopoints = [[-20,0.25], [-10,0.25], [-5ℯ/2,0.25], [-2ℯ,0.25], [-3ℯ/2,0.25], [-5ℯ/4,0.25], [-9ℯ/8,0.25], 
[-ℯ,0.25], [-7ℯ/8,0.25], [-3ℯ/4,0.25], [-1,1-1/ℯ], [-1,0.75], [-1,7/8], [-1,15/16], [-1,0.5] ]

# ╔═╡ d8613d71-af5e-4eb9-b065-b6467eae28b0
renatopoints = [[-3ℯ/2,4.5], [-4.05,4.75], [-2ℯ,2], [-ℯ,0.0625], [-3ℯ/4,0.5], [-ℯ/2,1],
	[-ℯ/4,0.25], [-ℯ/4,0.75], [-ℯ/4,1], [-ℯ/4,1.75], [-ℯ/8,2], [-ℯ/2,1.5], [-3ℯ/4,4], [-5ℯ/8,5], [-4ℯ/3,5.5], [-9ℯ/8,4.75], [-3ℯ/16,5.5]]

# ╔═╡ ee2b8e91-d405-4045-9287-e72daf20c7fd
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	cm = :grays #cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
			
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.5:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025), 
		maxiterations=100, colormap=cm
	)

	lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)
	scatter!(ax, Point2f.(marcopoints), color=1:length(marcopoints), colormap=:prism)
	scatter!(ax, Point2f.(renatopoints), color=:green, marker=:cross, markersize=12)

	fig
end

# ╔═╡ 7a0d58fd-7238-4f0c-accf-9e429ca4d325
md"""
## Software Julia packages
"""

# ╔═╡ a07e4983-d801-47ba-9531-16bdfb1f5f26
const Gr = SDDGraphics

# ╔═╡ cd230555-d920-4eb7-863c-e4ef586e84a0
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	cm = Gr.reverse(:grays) #cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
			
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.5:ymax, ygridcolor=gridc)
	
	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025) || real(z)<-16, 
		maxiterations=60, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)

	fs = 0.2 # FontSize
	ms = 16 # MarkerSize

	mypoints = [[-4ℯ/3,5.5], [-5ℯ/8,5], [-9ℯ/8,4.75]] # [-3ℯ/16,5.5]
	scatter!(ax, Point2f.(mypoints), color=:red, marker=:cross, markersize=ms)
	text!(ax, Point2f(-4ℯ/3,5.5), text="(-4ℯ/3,5.5)",fontsize=fs, color=:red, 
		  align = (:center,:bottom), markerspace=:data, offset=(0.0,0.05))
	text!(ax, Point2f(-9ℯ/8,4.75), text="(-9ℯ/8,4.75)",fontsize=fs, color=:red, 
		  align = (:center,:top), markerspace=:data, offset=(0.0,-0.05))
	text!(ax, Point2f(-5ℯ/8,5), text="(-5ℯ/8,5)",fontsize=fs, color=:red, 
		  align = (:center,:bottom), markerspace=:data, offset=(0.0,0.05))

	#save("Conjecture_plot.jpg", fig)
	
	fig
end

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
	a = 1+4im

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

	m = abs(g(a,a))
	#c1 = compassC(z->g(a,z), real(a) >= 0 ? a+1 : g(a,a)/m, m/2) # a+1.5a/abs(a)
	c1 = compassC(z->g(a,z), real(a) >= 0 ? a+2 : m, real(a) >= 0 ? 1 : 0.9m) # a+1.5a/abs(a)
	
	scatter!(ax2, [Point2f(real(a), imag(a))], color=:blue, markersize=4)
	scatter!(ax2, [Point2f(real(g(a,a)), imag(g(a,a)))], color=:green, markersize=8)
	scatter!(ax2, [Point2f(real(c1), imag(c1))], color=:red, markersize=8)
	
	fig
end

# ╔═╡ 8bcefc7f-d7ba-4566-a3be-7400f2797964
let
	a = 1

	Npix = 600
	xmin,xmax,ymin,ymax = -20,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 60
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,4Npix/3))
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

	c1 = compassC(z->g(a,z), real(a) >= 0 ? a+1 : g(a,a)/abs(g(a,a)), abs(g(a,a))/2)	
	
	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.75)], linewidth=0.5, markersize=4, iterations=100)
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

# ╔═╡ 70f58eb7-adc8-451c-94a6-11344bbb200b
let
	l,m = -20,0.25 # -4,4.5

	Npix = 400
	xmin,xmax,ymin,ymax = -12,2,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 64
	cm = cuherx
	
	fig = Figure(size=(3Npix/4,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))

	f = createflm(l,m)
	
	trappedpoints!(ax, f, xs, ys, hasescaped=z->real(z)<max(-64,min(-8,l/m*4)), maxiterations=maxits,
		colormap=cm #vermeerx #Gr.reverse(:vangogh)
	)

	c1 = findcritic1b(l, m, maxiterations=20)
	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	orbitpath!(ax, f, c1, colormap=[:red], linewidth=0.5, markersize=4, iterations=20)
	orbitpath!(ax, f, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c2), imag(c2))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c1), -imag(c1)), Point2f(real(c2), -imag(c2))], color=:green, markersize=8)	

	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)
	
	fig
end

# ╔═╡ aba0a028-b304-4a8d-91d3-f52693f34559
let
	l,m = -1,0.25 # -4,4.5

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 75
	cm = cuherx #vermeerx #Gr.reverse(:vangogh)
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,4Npix/3))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	f = createflm(l,m)
	f2 = createflm2(l,m)
	
	trappedpoints!(ax, f, xs, ys, maxiterations=maxits,
		hasescaped = z -> real(z) < -16, 
		colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#c1 = findcritic1b(l, m, maxiterations=20)
	c1 = findcriticBif(0.001,2.999, l, m, maxiterations=20)
	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)
	c3 = findcriticN(3, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	text!(ax, Point2f(real(c1),imag(c1)), text=L"c_1", color=:red, align=(:left,:center),
		  fontsize=0.5, markerspace=:data, offset=(0.25,0))
	text!(ax, Point2f(real(c2),imag(c2)), text=L"c_2", color=:blue, align=(:right,:center),
		  fontsize=0.5, markerspace=:data, offset=(-0.25,0))
	text!(ax, Point2f(real(c3),imag(c3)), text=L"c_3", color=:green, align=(:right,:center),
		  fontsize=0.5, markerspace=:data, offset=(-0.25,0))

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.75)], linewidth=0.5, markersize=4, iterations=100)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m025.jpg", fig)
	
	fig
end

# ╔═╡ ab535ab3-bb51-496d-b984-aef9a1b991d1
let
	Npix = 1200
	xmin,xmax,ymin,ymax = -2,4.5,-1.25,1.25
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
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

	save("Sienra_Leaf_0.jpg", fig)
	
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
		#seed = a -> compassC(z->g(a,z), real(a) >= 0 ? a+1 : g(a,a)/abs(g(a,a)), abs(g(a,a))/2),
		seed = a -> compassC(z->g(a,z), real(a) >= 0 ? a+2 : abs(g(a,a)), real(a) >= 0 ? 1 : 0.9abs(g(a,a))),
		#hasescaped = (a,z) -> real(z)>36,
		#hasescaped = (a,z) -> abs2(z-f(a,z))<0.001,
		hasescaped = (a,z) -> real(z)>64 || abs2(z-f(a,z))<0.001,
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

# ╔═╡ 99dc3a03-9aa0-4f52-8ed7-ffb4d4e5b538
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
	cm = cuherx #Gr.reverse(:cubehelix) # vermeerx #Gr.reverse(:vangogh)

	img = imgmandelbrot(
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcritic1b(real(lm), imag(lm), maxiterations=20, ε=0.001),
		hasescaped=(c,z)->abs(real(z))>16,
		maxiterations=maxits, colormap=cm
	)

	#save("Mandelbrot1_hires.jpg", img)
	
	fig = Figure(size=(Npix+200,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), backgroundcolor = RGBA(1,1,1,0),
		xtickformat = values -> ["$((xmax-xmin)*(value/Npix)+xmin)" for value in values], xticks=0:Npix/6:Npix, xgridcolor=RGBA(0.2,0.2,0.2,0.2),
		ytickformat = values -> ["$((ymax-ymin)*(value/Npix)+ymin)" for value in values], 
		yticks=0:Npix/6:Npix, ygridcolor=RGBA(0.2,0.2,0.2,0.2)) 
	#, yticks=ymin:1:ymax ) #, limits=(xmin,xmax,ymin,ymax), )

	#=mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcritic1b(real(lm), imag(lm), maxiterations=20, ε=0.001),
		hasescaped=(c,z)->real(z)<-16, # ||real(z)>16, 
		maxiterations=maxits, colormap=cm
	)=#
	image!(ax, rotr90(img), interpolate=true, fxaa=true, ssao=true, depth_shift=1)

	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)
	
	#save("Madelbrot1_plot.jpg", fig)
	
	fig
end

# ╔═╡ 6b8acca5-4fd3-4812-a6c9-7f75f6b17b64
mandelbrot(
	(c,z)->z^2+c, -2.01:0.002:0.51,-1.16:0.002:1.16, seed=0, colormap=cuherx, 
	hasescaped = (c,z) -> abs2(z)>4, maxiterations=48,
	depth_shift=1, interpolate=true, fxaa=true, ssao=true, # transparency=false, alpha=0.9, 
	axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16), xticks=-2:0.25:0.5, yticks=-1:0.25:1, xgridcolor=RGBA(0.25,0.25,0.25,0.25), ygridcolor=RGBA(0.25,0.25,0.25,0.25), backgroundcolor=RGBA(1,1,1,0)))

# ╔═╡ 35911590-c861-465d-921c-de93b66de7e7
let
	Npix = 500
	xmin,xmax,ymin,ymax = -21,0,0,7
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuherx #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(2Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> abs2(z-f2(real(lm), imag(lm), z))<0.0000025 && real(f(real(lm), imag(lm), z))<-16,
				#stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025) && real(f(real(lm), imag(lm), z))<-16, 
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)
	scatter!(ax, [Point2f(-20,0.25)], markersize=12, color=:red)
	text!(Point2f(-20,0.25), text="(-20,0.25)", color=:red, align = (:left,:bottom), fontsize=0.5,
		markerspace= :data)	

	Colorbar(fig[2, 1], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = false, flipaxis = false)

	#save("Bdomain_plot_rev.jpg", fig)
	
	fig
end

# ╔═╡ 0947a4b6-0309-4fea-bcab-7a3f0d036974
cuhex = :cubehelix #pushfirst!(push!(RGBA.(colorschemes[:cubehelix].colors), RGBA(1,1,1,0.1)),RGBA(0,0,0,0.9))

# ╔═╡ 74aa3553-400b-4df5-b053-7ee27338d8b6
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuhex #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025), 
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	fs = 0.096
	ms = 12
	
	lines!(ax,[-ℯ,0],[0,1],color=:cadetblue,linewidth=1.75)
	text!(Point2f(-0.18,0.74), text=L"\mathcal{L}_{-1}",
		  color=:dodgerblue4, align = (:center,:bottom), fontsize=fs, markerspace= :data, offset=(0.025,0.025))	
	
	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    	label = "Iterations", vertical = true, flipaxis = true)	
	
	#save("BulbsPer2_plot.jpg", fig)
	
	fig
end

# ╔═╡ 2176a612-f784-441b-8d44-ba132d5883b8
let
	Npix = 500
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuhex #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f2(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f2(real(lm), imag(lm), z), z, 2, ε=0.0000025), 
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	fs = 0.096
	ms = 12
	
	lines!(ax,[-ℯ,0],[0,1],color=:cadetblue,linewidth=1.75)
	text!(Point2f(-0.18,0.74), text=L"\mathcal{L}_{-1}",
		  color=:dodgerblue4, align = (:center,:bottom), fontsize=fs, markerspace= :data, offset=(0.025,0.025))	
	
	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    	label = "Iterations", vertical = true, flipaxis = true)	
	
	#save("BulbsPer2Finite_plot.jpg", fig)
	
	fig
end

# ╔═╡ a46480b9-1b34-4c81-82e1-be4af1ed3e05
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuhex #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025) && abs2(f(real(lm), imag(lm), z))<36, 
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	fs = 0.096
	ms = 12
	
	lines!(ax,[-ℯ,0],[0,1],color=:cadetblue,linewidth=1.75)
	text!(Point2f(-0.18,0.74), text=L"\mathcal{L}_{-1}",
		  color=:dodgerblue4, align = (:center,:bottom), fontsize=fs, markerspace= :data, offset=(0.025,0.025))	
	
	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    	label = "Iterations", vertical = true, flipaxis = true)	
	
	#save("BulbsPer2Finite_plot.jpg", fig)
	
	fig
end

# ╔═╡ 08e7bb0a-e954-46e0-9be8-eec3fce7d864
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 100
	cm = cuhex #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> real(f(real(lm), imag(lm), z))<-16,
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	fs = 0.096
	ms = 12
	
	lines!(ax,[-ℯ,0],[0,1],color=:cadetblue,linewidth=1.75)
	text!(Point2f(-0.18,0.74), text=L"\mathcal{L}_{-1}",
		  color=:dodgerblue4, align = (:center,:bottom), fontsize=fs, markerspace= :data, offset=(0.025,0.025))	
	
	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    	label = "Iterations", vertical = true, flipaxis = true)	
	
	#save("BulbBaker_plot.jpg", fig)
	
	fig
end

# ╔═╡ 3d076c9b-263b-4d62-a17e-afbd7fee2f21
let
	Npix = 400
	xmin,xmax,ymin,ymax = -3.5,0,0,2
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuhex #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,5Npix/8))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:0.5:xmax, xgridcolor=gridc,
			yticks=ymin:0.25:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025) && abs2(f(real(lm), imag(lm), z))<36, 
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	fs = 0.096 # FontSize
	ms = 12 # MarkerSize
	
	lines!(ax,[-ℯ,0],[0,1],color=:cadetblue,linewidth=1.75)
	text!(Point2f(-0.18,0.74), text=L"\mathcal{L}_{-1}",
		  color=:dodgerblue4, align = (:center,:bottom), fontsize=fs, markerspace= :data, offset=(0.025,0.025))	
	
	bifs1 = [Point2f(-1,0), Point2f(-1,1-1/ℯ), Point2f(-1,0.909075), Point2f(-1,1), Point2f(-1,1.9)]
	scatterlines!(ax, bifs1, markersize=ms, color=:gold)
	text!(ax, Point2f(-1,0), text="(-1,0)", fontsize=fs, color=:gold, 
		  align = (:left,:bottom), markerspace=:data, offset=(0.025,0.025))
	scatter!(ax, Point2f(-1,1-1/ℯ), markersize=ms, color=:black)
	text!(ax, Point2f(-1,1-1/ℯ), text="(-1,1-1/ℯ)", fontsize=fs, color=:black, 
		  align = (:left,:center), markerspace=:data, offset=(0.025,0.0))
	text!(ax, Point2f(-1,0.90975), text="(-1,0.090975)", fontsize=fs, color=:gold, 
		  align = (:left,:center), markerspace=:data, offset=(0.025,0.0))
	scatter!(ax, Point2f(-1,1), markersize=ms, color=:white)
	text!(ax, Point2f(-1,1), text="(-1,1)", fontsize=fs, color=:white, 
		  align = (:right,:center), markerspace=:data, offset=(-0.025,0.0))
	text!(ax, Point2f(-1,1.9), text="(-1,1.9)",fontsize=fs, color=:gold, 
		  align = (:left,:center), markerspace=:data, offset=(0.025,0.0))
	
	bifs2 = [Point2f(-1,0.25), Point2f(-3ℯ/4,0.25), Point2f(-2.079445,0.25), Point2f(-3.3,0.25)]
	scatterlines!(ax, bifs2, markersize=ms, color=:red)
	text!(Point2f(-1,0.25), text="(-1,1/4)", fontsize=fs, color=:red, 
		  align = (:left,:center), markerspace=:data, offset=(0.025,0.0))			
	text!(ax, Point2f(-3ℯ/4,0.25), text="(-3ℯ/4,1/4)", fontsize=fs, color=:red, 
		  align = (:left,:bottom), markerspace=:data, offset=(0.025,0.025))			
	text!(Point2f(-2.079445,0.25), text="(-2.079445,1/4)", fontsize=fs, color=:red, 
		  align = (:right,:top), markerspace=:data, offset=(-0.025,-0.025))	
	text!(Point2f(-3.3,0.25), text="(-3.3,1/4)", fontsize=fs, color=:red, 
		  align = (:left,:bottom), markerspace=:data, offset=(-0.025,0.025))			
	
	Colorbar(fig[2, 1], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = false, flipaxis = false)
	
	#save("Bifurcations_plot.jpg", fig)
	
	fig
end

# ╔═╡ 7e516588-5ad2-498b-89a6-91a376600ab2
cuherx_dark = pushfirst!(Gr.reverse(:cubehelix),RGB(0,0,0))

# ╔═╡ Cell order:
# ╟─afe66166-e9ca-11f0-bad5-59bd6959ad58
# ╠═f6b33556-cda7-4e17-8b1d-22801fb00c85
# ╠═b4c616d2-e720-465d-927d-e8534dc6cbec
# ╠═8e255a87-fe1c-4ab0-9330-1128cc317744
# ╠═e9a1d66a-abe1-4451-96ce-eccc6e887731
# ╠═c7fa0c18-f71a-4afa-8541-5f00dbfad2b5
# ╠═7d58b085-2f4d-4f2f-a392-5074b98d7663
# ╠═dcc54357-7342-4978-a0e3-ae25b3161158
# ╠═04d2ed2f-de6c-4e5a-a5de-c848b6de8a87
# ╠═729aaaf9-28f3-4274-9a26-12fc6f80f074
# ╠═4662f2ab-8b41-4211-98e0-58778d61d007
# ╠═f09e37b5-b0af-46c3-bf25-51017bf48ffa
# ╠═a2cf3c8b-c5c8-4df5-951e-62a7954e1083
# ╠═9bceda1a-3a79-4e35-95b2-2572d5577993
# ╠═cada7f65-ccb1-4ec3-98e6-2370aac392cc
# ╟─7087d5f4-86a4-4ed6-9970-3361ab584c79
# ╠═559b6103-c4f9-4edb-9f30-fc51ca2184cf
# ╠═8a8619e8-0236-4c37-bd09-65d82110fc2a
# ╠═502034db-8faa-41ae-92be-c150305a1105
# ╠═7e62773c-550b-483e-bf86-efed99abbbc6
# ╠═3031e81a-5cd4-47f5-8d24-0d407ea29327
# ╟─2f247a63-5943-49e7-a71e-f6d94af13181
# ╠═75927946-3b6a-4a06-ba52-669066a434f5
# ╠═6dce489a-ec3d-4907-b2b6-01e9043fe394
# ╠═1e61565e-072a-4d53-a573-41f91f191074
# ╟─fdbde74c-6ede-422c-955f-de5e02ec8b22
# ╠═978cc14e-3383-450b-a46f-0376baa84b8c
# ╠═781de362-0f5f-44e5-87cc-d77345617ca7
# ╠═d2a17754-e606-48e9-82a1-a350091538d2
# ╠═0bce24f6-9278-464e-9d0c-0af471089dec
# ╠═e3130dbc-4b0d-4fa0-80a2-2980fd875d0b
# ╠═e8e876c1-a196-406a-b444-a8618a208690
# ╠═639e9512-7281-4bc7-8bfc-76f4409306c6
# ╠═6b7870f5-a79d-416c-9a2d-6bd8a38e5fc5
# ╠═4cc91117-7c7b-4e70-820c-aff85860ab02
# ╠═f58338d6-a72a-4e23-8814-d3a88a5bb73f
# ╠═b02bddf3-40c1-4cb9-bdee-7de17456bba1
# ╟─208f363e-dc10-4307-aa27-3299ba2057a5
# ╠═8bcefc7f-d7ba-4566-a3be-7400f2797964
# ╟─e410e775-a417-40a2-8dbf-4c9c5050eef7
# ╟─70f58eb7-adc8-451c-94a6-11344bbb200b
# ╟─7273ed9e-5dfb-475b-b2a1-a65bfb64ab87
# ╠═aba0a028-b304-4a8d-91d3-f52693f34559
# ╠═e0ff84af-11fd-4c81-8792-8c45e20984b2
# ╟─ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
# ╠═ab535ab3-bb51-496d-b984-aef9a1b991d1
# ╠═d5872959-a8a4-443a-909a-265eddb9b5a1
# ╠═cdf5ac1b-df7f-4270-8c67-1340e20861b4
# ╠═99dc3a03-9aa0-4f52-8ed7-ffb4d4e5b538
# ╟─0fd2c626-52e9-4f8c-8eac-64a2047d9729
# ╟─6b8acca5-4fd3-4812-a6c9-7f75f6b17b64
# ╟─d91916f8-37a6-4494-ac61-6df7aba6a0c1
# ╟─6e0bf940-1d94-4efa-a0a0-6558ba5d67fe
# ╠═35911590-c861-465d-921c-de93b66de7e7
# ╟─65365633-6555-4f3d-945b-66275dcf8c7b
# ╠═ffbacbc1-1b56-4b1f-93f1-0dce2f44bc9f
# ╠═74aa3553-400b-4df5-b053-7ee27338d8b6
# ╟─249dda01-d980-4cc6-86f0-1df2a31dd7d6
# ╟─07765e47-4d63-4df8-ab8f-906d71a75338
# ╟─2b61308c-3efa-47a9-88d8-0ea870ca3d77
# ╠═2176a612-f784-441b-8d44-ba132d5883b8
# ╠═855b2724-3129-4166-b170-99f6ae6cbc37
# ╠═a46480b9-1b34-4c81-82e1-be4af1ed3e05
# ╟─ab526193-7dfd-4e82-a397-beb7e2fd6c93
# ╟─e23f79d7-5899-409b-b721-4bb06aac37b4
# ╠═08e7bb0a-e954-46e0-9be8-eec3fce7d864
# ╟─40f1cf47-2b5c-4d9e-bc8b-10ca781dffcc
# ╠═3d076c9b-263b-4d62-a17e-afbd7fee2f21
# ╟─05efe4d3-4126-449a-ac48-1799c4c2f28a
# ╟─de1c003f-d951-47cf-acc3-beaea021298b
# ╟─2c9f1a59-ad76-4476-8d32-7e158dcd2541
# ╠═ee95a23b-177b-474d-a517-5cc9204831dd
# ╠═fab68d9a-311d-451f-9289-caaec83569dc
# ╠═d8613d71-af5e-4eb9-b065-b6467eae28b0
# ╠═ee2b8e91-d405-4045-9287-e72daf20c7fd
# ╠═cd230555-d920-4eb7-863c-e4ef586e84a0
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
