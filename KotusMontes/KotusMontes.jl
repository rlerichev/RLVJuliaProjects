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
using Colors, ColorSchemes, Images

# ╔═╡ df10e677-8b09-473e-8b15-3ee063496f34
using Makie, GLMakie

# ╔═╡ 4294274f-8874-4840-a241-0cce4969d5dd
using SDDCore, SDDGraphics, SDDGeometry, SDD

# ╔═╡ d52d0e77-9a16-4163-b2ef-10779cd1cdff
include("../common/src/nbutils.jl")

# ╔═╡ afe66166-e9ca-11f0-bad5-59bd6959ad58
md"""
# Dynamics of $f_{\lambda,\mu}$

### Function and derivatives

$f_{\lambda,\mu}(z)=\lambda e^z + \frac{\mu}{z}$

$f'_{\lambda,\mu}(z)=\lambda e^z - \frac{\mu}{z^2}$

$f''_{\lambda,\mu}(z)=\lambda e^z + 2\frac{\mu}{z^3}$

$f^2_{\lambda,\mu}(z)=\lambda e^{\lambda e^z + \frac{\mu}{z}} + \frac{\mu}{\lambda e^z + \frac{\mu}{z}}=\lambda e^{\lambda e^z + \frac{\mu}{z}} + \frac{\mu z}{\lambda z e^z + \mu}$
"""

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
## Critical points: Intersection of curves

 $\lambda,\mu\in\mathbb{R}$, $\lambda<0$ & $\mu>0$.

$f'(z)=0 \implies z^2 e^z = \frac{\mu}{\lambda}$

1.

$\implies |z^2 e^z| = (x^2+y^2) e^x = \big|\frac{\mu}{\lambda}\big|$

$\implies y = y(x)= \pm\sqrt{\big|\frac{\mu}{\lambda}\big|\frac{1}{e^x}-x^2}$

2.

$\implies \arg(z^2 e^z) = 2\arg(x+yi)+y+2\pi k = \arg\big(\frac{\mu}{\lambda}\big)=-\pi$

$\implies \tan^{-1}\big(\frac{y}{x}\big) = \frac{-\pi - y -2\pi k}{2}$

$\implies x = x(y) = \frac{y}{\tan\big(\frac{-\pi - y -2\pi k}{2}\big)} = -\frac{y}{\tan\big(\frac{y+\pi}{2}\big)}$

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

# ╔═╡ dbbdd5ab-3e34-4419-b4d4-a6d2d9d2f558
let
	l,m = -20,0.25
	xmin,xmax,ymin,ymax = -0.8,0.8,0,1.8
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))
	#vlines!(ax, [-0.5m/l], color=:black, linestyle=:dash)
	xs = xmin:0.001:xmax
	lines!(ax, xs, cy.(xs,l,m), color=:cadetblue, linewidth=4)
	ys = ymin:0.001:ymax
	lines!(ax, cx.(ys,l,m), ys, color=:orange, linewidth=4)
	fig
end

# ╔═╡ 63987b6c-b0d6-4353-98af-c769dd7ca751
let
	l,m = -20,0.25
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(-12,2,0,12))
	xs = -12:0.02:2
	lines!(ax, xs, cy.(xs,l,m), color=:cadetblue, linewidth=4)
	ys = 0:0.02:12
	lines!(ax, cx.(ys,l,m), ys, color=:orange, linewidth=4)
	hlines!(ax, [pi,2pi,3pi], color=:black, linestyle=:dash)
	fig
end

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

# ╔═╡ 208f363e-dc10-4307-aa27-3299ba2057a5
md"""
## Dynamic plane
"""

# ╔═╡ 369f53e9-4ea9-4462-8297-ebb3fc381ffb
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

# ╔═╡ edf8fc8f-b6bd-4f29-af08-00039b30acd1
let
	l,m = -20,0.25 #-1.5,5 #-2,0.5 #0.1,4 #-4,4.5 #-4,1 #-3,0.1 #-0.25,1 # -1,0.25
	
	xmin,xmax,ymin,ymax = -12,2,-6,6
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))

	f = createflm(l,m)
	
	trappedpoints!(ax, f, xs, ys, hasescaped=z->stops(f,z,2) || real(z)<-12, maxiterations=20, colormap=:vangogh)

	c1 = findcritic1b(l, m, maxiterations=20)
	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	orbitpath!(ax, f, c1, colormap=[:red], linewidth=0.5, markersize=4, iterations=20)
	orbitpath!(ax, f, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c2), imag(c2))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c1), -imag(c1)), Point2f(real(c2), -imag(c2))], color=:green, markersize=8)	
	
	fig
end

# ╔═╡ ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
md"""
## Parameter plane
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
let
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
end

# ╔═╡ de1c003f-d951-47cf-acc3-beaea021298b
let
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
end

# ╔═╡ 7a0d58fd-7238-4f0c-accf-9e429ca4d325
md"""
## Software Julia packages
"""

# ╔═╡ a07e4983-d801-47ba-9531-16bdfb1f5f26
const Gr = SDDGraphics

# ╔═╡ 99dc3a03-9aa0-4f52-8ed7-ffb4d4e5b538
let
	xmin,xmax,ymin,ymax = -10,0,0,10
	Δx,Δy = (xmax-xmin)/920, (ymax-ymin)/920
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(800,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:1:xmax, yticks=ymin:1:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcritic1b(real(lm), imag(lm), maxiterations=20, ε=0.001),
		hasescaped=(c,z)->real(z)<-16, # ||real(z)>16, 
		maxiterations=80, colormap=Gr.reverse(:vangogh)
	)

	save("test.jpg", fig)
	
	fig
end

# ╔═╡ 6b8acca5-4fd3-4812-a6c9-7f75f6b17b64
mandelbrot((c,z)->z^2+c, -2.01:0.005:0.51,-1.16:0.005:1.16, seed=0, colormap=Gr.reverse(:vangogh),
	hasescaped = (c,z) -> abs2(z)>4, maxiterations=60, axis=(;aspect=DataAspect()))

# ╔═╡ c14eacb3-f173-4609-b750-413027a296ca
prismx = push!(RGBA.(deepcopy(colorschemes[:prism].colors)), RGBA(1,1,1,0.25))

# ╔═╡ fd700d57-b159-4a0e-87e6-6b7361254b60
vermeerx = push!(RGBA.(Gr.reverse(:vermeer)), RGBA(0,0,0,0.5))

# ╔═╡ 70f58eb7-adc8-451c-94a6-11344bbb200b
let
	l,m = -20,0.25
	
	xmin,xmax,ymin,ymax = -12,2,-6,6
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))

	f = createflm(l,m)
	
	trappedpoints!(ax, f, xs, ys, hasescaped=z->real(z)<-12, maxiterations=120,
		colormap=vermeerx #Gr.reverse(:vangogh)
	)

	c1 = findcritic1b(l, m, maxiterations=20)
	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	orbitpath!(ax, f, c1, colormap=[:red], linewidth=0.5, markersize=4, iterations=20)
	orbitpath!(ax, f, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c2), imag(c2))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c1), -imag(c1)), Point2f(real(c2), -imag(c2))], color=:green, markersize=8)	
	
	fig
end

# ╔═╡ 2b61308c-3efa-47a9-88d8-0ea870ca3d77
mandelbrot((c,z)->z^2+c, -2.01:0.005:0.51,-1.16:0.005:1.16, seed=0, colormap=vermeerx, #Gr.reverse(:vermeer),
	hasescaped = (c,z) -> stops(z->z^2+c, z, 2, ε=0.00001) || abs2(z)>4, maxiterations=240, axis=(;aspect=DataAspect()))

# ╔═╡ 2382c3cc-266c-4f82-9c2f-daf60ded2ee9
let
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
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.00001), 
		maxiterations=40, colormap=vermeerx #:vangogh
	)

	fig
end

# ╔═╡ 74ebf7f1-da95-40fe-bc67-2e2dee196df6
let
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:1:xmax, yticks=ymin:1:ymax)

	mandelbrot!(ax,
		(lm,z) -> f2(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> abs(real(f(real(lm), imag(lm), z)))>16,
		#hasescaped = (lm,z) -> stops(z->f2(real(lm), imag(lm), z), z, 2, ε=0.00001),
		#|| real(f(real(lm), imag(lm), z))<-16, 
		maxiterations=40, colormap=vermeerx #:vangogh
	)

	fig
end

# ╔═╡ 05efe4d3-4126-449a-ac48-1799c4c2f28a
let
	xmin,xmax,ymin,ymax = -5,-3,4,5
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:1:xmax, yticks=ymin:1:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.00001), # || abs2(z)>100, 
		maxiterations=320, colormap=vermeerx #:vangogh
	)

	fig
end

# ╔═╡ 2c9f1a59-ad76-4476-8d32-7e158dcd2541
let
	xmin,xmax,ymin,ymax = -5,0,0,2
	Δx,Δy = (xmax-xmin)/800, (ymax-ymin)/800
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1000,600))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax), xticks=xmin:1:xmax, yticks=ymin:1:ymax)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.000001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.00001), # || real(z)<-16, 
		maxiterations=120, colormap=vermeerx #Gr.reverse(:vermeer) #prismx #:vangogh
	)

	#lines!(ax, [-ℯ,0], [0,2/ℯ], color=:red)

	fig
end

# ╔═╡ Cell order:
# ╟─afe66166-e9ca-11f0-bad5-59bd6959ad58
# ╠═dcc54357-7342-4978-a0e3-ae25b3161158
# ╠═04d2ed2f-de6c-4e5a-a5de-c848b6de8a87
# ╠═729aaaf9-28f3-4274-9a26-12fc6f80f074
# ╠═4662f2ab-8b41-4211-98e0-58778d61d007
# ╠═f09e37b5-b0af-46c3-bf25-51017bf48ffa
# ╠═9bceda1a-3a79-4e35-95b2-2572d5577993
# ╠═cada7f65-ccb1-4ec3-98e6-2370aac392cc
# ╟─8a8619e8-0236-4c37-bd09-65d82110fc2a
# ╠═502034db-8faa-41ae-92be-c150305a1105
# ╠═7e62773c-550b-483e-bf86-efed99abbbc6
# ╠═3031e81a-5cd4-47f5-8d24-0d407ea29327
# ╠═dbbdd5ab-3e34-4419-b4d4-a6d2d9d2f558
# ╠═63987b6c-b0d6-4353-98af-c769dd7ca751
# ╟─2f247a63-5943-49e7-a71e-f6d94af13181
# ╠═75927946-3b6a-4a06-ba52-669066a434f5
# ╠═6dce489a-ec3d-4907-b2b6-01e9043fe394
# ╠═1e61565e-072a-4d53-a573-41f91f191074
# ╠═fdbde74c-6ede-422c-955f-de5e02ec8b22
# ╠═978cc14e-3383-450b-a46f-0376baa84b8c
# ╠═781de362-0f5f-44e5-87cc-d77345617ca7
# ╠═d2a17754-e606-48e9-82a1-a350091538d2
# ╠═0bce24f6-9278-464e-9d0c-0af471089dec
# ╠═e3130dbc-4b0d-4fa0-80a2-2980fd875d0b
# ╠═e8e876c1-a196-406a-b444-a8618a208690
# ╠═6b7870f5-a79d-416c-9a2d-6bd8a38e5fc5
# ╠═4cc91117-7c7b-4e70-820c-aff85860ab02
# ╟─208f363e-dc10-4307-aa27-3299ba2057a5
# ╠═70f58eb7-adc8-451c-94a6-11344bbb200b
# ╠═369f53e9-4ea9-4462-8297-ebb3fc381ffb
# ╠═edf8fc8f-b6bd-4f29-af08-00039b30acd1
# ╟─ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
# ╠═99dc3a03-9aa0-4f52-8ed7-ffb4d4e5b538
# ╠═6b8acca5-4fd3-4812-a6c9-7f75f6b17b64
# ╠═d91916f8-37a6-4494-ac61-6df7aba6a0c1
# ╠═2b61308c-3efa-47a9-88d8-0ea870ca3d77
# ╠═6e0bf940-1d94-4efa-a0a0-6558ba5d67fe
# ╠═2382c3cc-266c-4f82-9c2f-daf60ded2ee9
# ╠═74ebf7f1-da95-40fe-bc67-2e2dee196df6
# ╠═05efe4d3-4126-449a-ac48-1799c4c2f28a
# ╠═de1c003f-d951-47cf-acc3-beaea021298b
# ╠═2c9f1a59-ad76-4476-8d32-7e158dcd2541
# ╠═d52d0e77-9a16-4163-b2ef-10779cd1cdff
# ╠═7a0d58fd-7238-4f0c-accf-9e429ca4d325
# ╠═7af6c311-4ff3-4f71-84f0-8d9ef2d591c8
# ╠═df10e677-8b09-473e-8b15-3ee063496f34
# ╠═4b8a3a21-889f-4189-a6c5-e06014e1e8e0
# ╠═4294274f-8874-4840-a241-0cce4969d5dd
# ╠═a07e4983-d801-47ba-9531-16bdfb1f5f26
# ╠═c14eacb3-f173-4609-b750-413027a296ca
# ╠═fd700d57-b159-4a0e-87e6-6b7361254b60
