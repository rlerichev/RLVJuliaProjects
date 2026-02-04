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
using Makie, GLMakie, CairoMakie

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

# ╔═╡ 72035530-a6b7-4c0e-8420-a4b2617d739e
md"""
## Attractor

$f_{\lambda,\mu}(z)=\lambda e^z+\frac{\mu}{z}=z$

$|f'_{\lambda,\mu}(z)|=|\lambda e^z-\frac{\mu}{z^2}|<1.$

$\implies$

$\lambda z e^z+\mu=z^2$

$|\lambda z^2 e^z-\mu|<|z^2|.$

$\implies$

$|\lambda z^2 e^z-\mu|^2<|\lambda z e^z+\mu|^2.$

$\implies$

$(\lambda z^2 e^z-\mu)(\overline{\lambda z^2 e^z-\mu})<(\lambda z e^z+\mu)(\overline{\lambda z e^z+\mu}).$

$\implies$

$(\lambda z^2 e^z-\mu)(\overline{\lambda z^2 e^z-\mu})<(\lambda z e^z+\mu)(\overline{\lambda z e^z+\mu}).$

 $z\neq 0$, since $0$ is the pole, then

$\frac{\lambda e^z}{z}+\frac{\mu}{z^2}=1$
"""

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
#=let
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
end=#

# ╔═╡ 63987b6c-b0d6-4353-98af-c769dd7ca751
#=let
	l,m = -20,0.25
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(-12,2,0,12))
	xs = -12:0.02:2
	lines!(ax, xs, cy.(xs,l,m), color=:cadetblue, linewidth=4)
	ys = 0:0.02:12
	lines!(ax, cx.(ys,l,m), ys, color=:orange, linewidth=4)
	hlines!(ax, [pi,2pi,3pi], color=:black, linestyle=:dash)
	fig
end=#

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

# ╔═╡ 12e25b26-5f8b-4bd8-9269-9b4672c66bf8
md"""
Bifurcation parameter

$(\lambda,\mu)=(-1,1-1/e)\approx(-1,0.632120558828)$

Parabolic fixed.
"""

# ╔═╡ 0e4944e6-0626-4644-bcc2-6c3f1a5c5b37
1-1/ℯ

# ╔═╡ a8ede40d-b3b4-47f7-9a3c-2a787dc771a6
f(-1,1-1/ℯ,-1), f´(-1,1-1/ℯ,-1)

# ╔═╡ 9a344400-1172-4666-b4cd-1daec8f7339f
md"""
$(\lambda,\mu)=(-1,0.771)$

Attracting period 2.
"""

# ╔═╡ 3e83fb22-774d-4c41-85c6-21567d542ecc
(1-1/ℯ+0.90975)/2

# ╔═╡ 01e1f581-df35-4179-9af2-3ce2e7a8ee20
md"""
Bifurcation parameter

$(\lambda,\mu)=(-1,0.90975)$

Parabolic period 2.
"""

# ╔═╡ d188982b-663f-4405-916d-461ae5e24ee5
let
	z1 = -0.12 + 1.49im
	z2 = conj(z1)

	f(-1,0.90975,z1), f´(-1,0.90975,z1)*f´(-1,0.90975,z2)
end

# ╔═╡ 5452962c-8c9b-4e27-bdb5-608259f16f49
md"""
$(\lambda,\mu)=(-1,15/16)=(-1,0.9375)$

Attracting period 2.
"""

# ╔═╡ a983327e-8782-45d8-8459-a10dd7ad209b
15/16

# ╔═╡ e28fc392-446f-443e-97ed-fd7de84efd5f
md"""
"Classical case"

$(\lambda,\mu)=(-1,1)$

Attracting period 2.
"""

# ╔═╡ d1e50b30-e677-4fe8-8cb7-e65aac67441b
md"""
Bifurcation parameter???

$(\lambda,\mu)=(-1,1.9)$

???
"""

# ╔═╡ 263e7c61-b903-4cf7-9c8b-2872fb0c9989
let
	z1 = -0.2385 + 2.8539im
	z2 = 0.7062 - 0.85885im	
	z3 = conj(z1)
	z4 = conj(z2)
	f(-1, 1.9, z1), f(-1, 1.9, z2), f´(-1, 1.9, z1)*f´(-1, 1.9, z2)
end

# ╔═╡ 9e603559-2a3d-4d05-829a-99ef2dac32fb
md"""
Bifurcation parameter

$(\lambda,\mu)=(-3e/4,1/4)\approx(-2.03871137,0.25)$

Parabolic fixed.
"""

# ╔═╡ 276212fc-720a-441d-ad15-0a3cc7c8c5a4
-3ℯ/4

# ╔═╡ 60e1123d-f204-4eca-a522-ebbc7d38e263
md"""
$(\lambda,\mu)=(-e,1/16)\approx(-2.71828,0.0625)$

Attracting fixed.
"""

# ╔═╡ a4fe1407-067c-4192-a6eb-844423a7a279
-ℯ, 1/16

# ╔═╡ 7bb26926-347f-4a43-8fa7-f3cff2bddfae
md"""
Bifurcation parameter

$(\lambda,\mu)=(-2.0794454,1/4)$

Parabolic period 2.
"""

# ╔═╡ d6c23140-9dba-49c8-8a64-4186ccf1ca9e
md"""
$(\lambda,\mu)=(-e,1/4)\approx(-2.71828,0.25)$

Attracting period 2.
"""

# ╔═╡ 795898a7-fabb-418f-aa19-ef98afb490f9
-4ℯ/3

# ╔═╡ 9bf7b086-d2d2-430a-a891-e43932c5c79d
md"""
$(\lambda,\mu)=(-3.3,0.25)$

Bifurcation point???
"""

# ╔═╡ 7e8d748e-50ee-48a5-9081-d7f0e400647d
md"""
$(\lambda,\mu)=(-5e/8,5)\approx(-1.698936,5)$
"""

# ╔═╡ 9693d1da-7886-4d7a-906d-004196bd338f
-5ℯ/8

# ╔═╡ 2bb4db46-4d9c-4eb0-8ea9-99d95d2804fa
md"""
$(\lambda,\mu)=(-9e/8,4.75)\approx(-3.05867)$
"""

# ╔═╡ ac30ebe5-915b-4f2a-9634-912d4d7ad369
-9ℯ/8

# ╔═╡ e17ec7eb-b95e-4e93-894c-3dab3cf26e10
md"""
$(\lambda,\mu)=(-4e/3,5.5)\approx(-3.624375,5.5)$
"""

# ╔═╡ 9f02bbfc-d1c4-43d7-a940-a6b6a6b7305b
-4ℯ/3

# ╔═╡ 0dc154e3-b602-42be-9b56-3325f1961581
#=let
	l,m = -20,0.25
	
	xmin,xmax,ymin,ymax = -12,2,-6,6
	Δx,Δy = (xmax-xmin)/800, (ymax-ymin)/800
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(800,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))

	f2 = createflm2(l,m)
	
	trappedpoints!(ax, f2, xs, ys, hasescaped=z->real(z)<-100, maxiterations=32,
		colormap=vermeerx #Gr.reverse(:vangogh)
	)

	c1 = findcritic1b(l, m, maxiterations=20)
	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	orbitpath!(ax, f2, c1, colormap=[:red], linewidth=0.5, markersize=4, iterations=20)
	orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c2), imag(c2))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c1), -imag(c1)), Point2f(real(c2), -imag(c2))], color=:green, markersize=8)	
	
	fig
end=#

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
#=let
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
end=#

# ╔═╡ 70a90845-4d5f-44e0-997b-4ff969b56346
#=let
	l,m = -20,0.25 #-1.5,5 #-2,0.5 #0.1,4 #-4,4.5 #-4,1 #-3,0.1 #-0.25,1 # -1,0.25
	
	xmin,xmax,ymin,ymax = -12,2,-6,6
	Δx,Δy = (xmax-xmin)/720, (ymax-ymin)/720
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax
	
	fig = Figure(size=(1200,800))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax))

	f2 = createflm2(l,m)
	
	trappedpoints!(ax, f2, xs, ys, hasescaped=z->stops(f2,z,2), maxiterations=32, colormap=vermeerx)

	c1 = findcritic1b(l, m, maxiterations=20)
	c2 = findcriticN(2, l, m, maxiterations=100, maxiterationscompass=16, ε=0.00000001)

	orbitpath!(ax, f2, c1, colormap=[:red], linewidth=0.5, markersize=4, iterations=20)
	orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c2), imag(c2))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c1), -imag(c1)), Point2f(real(c2), -imag(c2))], color=:green, markersize=8)	
	
	fig
end=#

# ╔═╡ ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
md"""
## Parameter plane
Let $c_1=c_1(\lambda,\mu),\,c_2=c_2(\lambda,\mu)\in\mathbb{C}$ be the critical points of $f_{\lambda,\mu}$ such that $Im(c_1)\in(0,\pi)$ and $Im(c_2)\in(\pi,3\pi)$.
"""

# ╔═╡ cdf5ac1b-df7f-4270-8c67-1340e20861b4
md"""
In black, the "Mandelbrot set" of $f_{\lambda,\mu}$:

$\mathbb{M}_1=\{(\lambda,\mu)\in\mathbb{R}^2|\,\,o(c_1,f_{\lambda,\mu})\,\,\mathrm{is\,bounded}\}.$

In colors, using the escape time:

$\mathbb{R}^2-\mathbb{M}_1.$
"""

# ╔═╡ 99dc3a03-9aa0-4f52-8ed7-ffb4d4e5b538
#=let
	Npix = 400
	xmin,xmax,ymin,ymax = -16,0,0,32
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
	cm = cuherx #Gr.reverse(:cubehelix) # vermeerx #Gr.reverse(:vangogh)

	#=img = imgmandelbrot(
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcritic1b(real(lm), imag(lm), maxiterations=20, ε=0.001),
		hasescaped=(c,z)-> abs(real(z))>16, # abs2(z)>10000,
		maxiterations=maxits, colormap=cm
	)=#

	#save("Mandelbrot1_hires.jpg", img)
	
	fig = Figure(size=(Npix+200,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), backgroundcolor = RGBA(1,1,1,0),
		#xtickformat = values -> ["$((xmax-xmin)*(value/Npix)+xmin)" for value in values], 
		xticks=xmin:1:xmax, xgridcolor=RGBA(0.2,0.2,0.2,0.2),
		#ytickformat = values -> ["$((ymax-ymin)*(value/Npix)+ymin)" for value in values], 			
		yticks=ymin:2:ymax, ygridcolor=RGBA(0.2,0.2,0.2,0.2)) 
	#, yticks=ymin:1:ymax ) #, limits=(xmin,xmax,ymin,ymax), )

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcritic1b(real(lm), imag(lm), maxiterations=20, ε=0.001),
		hasescaped=(c,z)->real(z)<-16, # ||real(z)>16, 
		maxiterations=maxits, colormap=cm
	)
	#image!(ax, rotr90(img), interpolate=true, fxaa=true, ssao=true, depth_shift=1)

	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)
	
	#save("Madelbrot1_plot.jpg", fig)
	
	fig
end=#

# ╔═╡ 7786dcfb-1d6e-419b-802c-e812a37361eb


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

# ╔═╡ 4c252064-eb11-485b-9607-01aaf84242bb
#specialmandelbrot(f::Function, xs::AbstractVector, xs::AbstractVector)

# ╔═╡ 735528fa-e490-47ab-bc6c-bbba9b2fdac8
md"""
$\mathbb{B}=\{(\lambda,\mu)|\,\,F_0\cap Crit(f_{\lambda,\mu})\neq\emptyset\}$
"""

# ╔═╡ 0ec76c85-d021-4a0b-8872-882c03f1a0d3
md"""
In black, the "Mandelbrot set" of $f_{\lambda,\mu}$

$\mathbb{M}'_1=\{(\lambda,\mu)|\,\,o(c_1,f_{\lambda,\mu})\,\mathrm{is\,bounded,\,or\,} f^n_{\lambda,\mu}(c_1)\rightarrow+\infty\}$
"""

# ╔═╡ 65365633-6555-4f3d-945b-66275dcf8c7b
md"""
Reference: Mandelbrot set.
"""

# ╔═╡ fef962ec-4550-4a2f-bb28-5df9c7ced028
#=mandelbrot((c,z)->z^2+c, -2.01:0.002:0.51,-1.16:0.002:1.16, seed=0, colormap=cuherx, #vermeerx,
	hasescaped = (c,z) -> !stops(z->z^2+c, z, 2, ε=0.00001) && abs2(z)>144, maxiterations=48, #  || 
	depth_shift=1, interpolate=true, fxaa=true, ssao=true, # transparency=false, alpha=0.9, 
	axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16), xticks=-2:0.25:0.5, yticks=-1:0.25:1, xgridcolor=RGBA(0.25,0.25,0.25,0.25), ygridcolor=RGBA(0.25,0.25,0.25,0.25), backgroundcolor=RGBA(1,1,1,0)))
=#

# ╔═╡ ffbacbc1-1b56-4b1f-93f1-0dce2f44bc9f
md"""
Let $z_0=z_0(\lambda,\mu)\in\mathbb{C}$ be the convergence point of the sequence $f_{\lambda,\mu}^{2n}(c_1)$, such that $|(f_{\lambda,\mu}^{2n})'(z_0)|<1$, or $z_0=0$, or $z_0\rightarrow\infty$.

In colors (coloring with $n\mapsto|f_{\lambda,\mu}^n(c_1)-z_0|<\varepsilon$):

$\mathbb{M}_2=\{(\lambda,\mu)\in\mathbb{R}^2|\,\, f_{\lambda,\mu}^{2n}(c_1)\rightarrow z_0\}.$

In black:

$\mathbb{R}^2-\mathcal{M}_2.$

In light blue, the line of parameters where $f_{\lambda,\mu}$ has fixed parabolic points with multiplier $-1$:

$\mathcal{L}_{-1}=\{(\lambda,\mu)\in\mathbb{R}^2|\,\,\exists z:\,\,f_{\lambda,\mu}(z)=z,\,\,f'_{\lambda,\mu}(z)=-1\}.$
"""

# ╔═╡ 249dda01-d980-4cc6-86f0-1df2a31dd7d6
md"""
Reference: in colors the Mandelbrot set's bulbs of period $1$ and $2$, union its exterior:

$\{z\in\mathbb{C}|\,\,q_c^{2n}(0)\rightarrow z_c\,\,\mathrm{or}\,\,q_c^{2n}(0)\rightarrow\infty\},$

where $q_c(z)=z^2+c$ and $z_c$ is an attracting fixed point of $q^2_c$.
"""

# ╔═╡ 07765e47-4d63-4df8-ab8f-906d71a75338
mandelbrot((c,z)->z^2+c, -2.01:0.01:0.51,-1.16:0.01:1.16, seed=0, colormap=:cubehelix, #vermeerx,
	hasescaped = (c,z) -> stops(z->z^2+c, z, 2, ε=0.00001) || abs2(z)>144, maxiterations=128,
	depth_shift=1, interpolate=true, fxaa=true, ssao=true, # transparency=false, alpha=0.9, 
	axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16), xticks=-2:0.25:0.5, yticks=-1:0.25:1, xgridcolor=RGBA(0.25,0.25,0.25,0.25), ygridcolor=RGBA(0.25,0.25,0.25,0.25), backgroundcolor=RGBA(1,1,1,0)))

# ╔═╡ 2b61308c-3efa-47a9-88d8-0ea870ca3d77
#=mandelbrot((c,z)->z^2+c, -2.01:0.002:0.51,-1.16:0.002:1.16, seed=0, colormap=cuherx, #vermeerx,
	hasescaped = (c,z) -> stops(z->z^2+c, z, 2, ε=0.00001) || abs2(z)>144, maxiterations=128, axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16)))
=#

# ╔═╡ d2598cd0-8ede-4558-bc38-b5503bab5f17
md"""
## Bulb Period 1
"""

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
mandelbrot((c,z)->z^2+c, -2.01:0.01:0.51,-1.16:0.01:1.16, seed=0, colormap=:cubehelix, #vermeerx,
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
md"""
$\mathbb{E}_k=\{(\lambda,\mu)|\,\,f^n_{\lambda,\mu}(c_k)\rightarrow\{0,\infty\}\}$
"""

# ╔═╡ 5cc6c94c-fccf-4251-aaf3-0629f396457f
md"""
---
"""

# ╔═╡ fab68d9a-311d-451f-9289-caaec83569dc
marcopoints = [[-20,0.25], [-10,0.25], [-5ℯ/2,0.25], [-2ℯ,0.25], [-3ℯ/2,0.25], [-5ℯ/4,0.25], [-9ℯ/8,0.25], 
[-ℯ,0.25], [-7ℯ/8,0.25], [-3ℯ/4,0.25], [-1,1-1/ℯ], [-1,0.75], [-1,7/8], [-1,15/16], [-1,0.5] ]

# ╔═╡ d8613d71-af5e-4eb9-b065-b6467eae28b0
renatopoints = [[-3ℯ/2,4.5], [-4.05,4.75], [-2ℯ,2], [-ℯ,0.0625], [-3ℯ/4,0.5], [-ℯ/2,1],
	[-ℯ/4,0.25], [-ℯ/4,0.75], [-ℯ/4,1], [-ℯ/4,1.75], [-ℯ/8,2], [-ℯ/2,1.5], [-3ℯ/4,4], [-5ℯ/8,5], [-4ℯ/3,5.5], [-9ℯ/8,4.75], [-3ℯ/16,5.5]]

# ╔═╡ ee2b8e91-d405-4045-9287-e72daf20c7fd
let
	Npix = 300
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

	plt = mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025), 
		maxiterations=100, colormap=cm
	)
	translate!(plt, 0, 0, -1)

	lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)
	scatter!(ax, Point2f.(marcopoints), color=1:length(marcopoints), colormap=:prism)
	scatter!(ax, Point2f.(renatopoints), color=:green, marker=:cross, markersize=12)

	#save("test.svg", fig)
	
	fig
end

# ╔═╡ 7a0d58fd-7238-4f0c-accf-9e429ca4d325
md"""
## Software Julia packages
"""

# ╔═╡ 92366d0e-43c7-49c9-b78a-07ced091dcee
#CairoMakie.activate!(type = "png")

# ╔═╡ 36239fa1-7865-4b54-bd6a-025c80bceae8
GLMakie.activate!()

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

# ╔═╡ 028a0597-e08f-4424-b379-b7f604a0b79e
let
	l,m = -1,1-1/ℯ # -4,4.5

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.75)], linewidth=0.5, markersize=4, iterations=300)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	scatter!(ax, Point2f(-1,0), color=:yellow, markersize=8, marker=:cross)	
	
	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m1-1e.jpg", fig)
	
	fig
end

# ╔═╡ 406d4a9f-539a-4510-aa56-6cb1e0dcfcf1
let
	l,m = -1,0.771

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=2, iterations=100)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m0771.jpg", fig)
	
	fig
end

# ╔═╡ fcb2bc42-4789-491b-b966-455569671de2
let
	l,m = -1,0.90975

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=2, iterations=200)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	scatter!(ax, [Point2f(-0.12,1.49), Point2f(-0.12,-1.49)], color=:yellow, markersize=8, marker=:cross)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m090975.jpg", fig)
	
	fig
end

# ╔═╡ 876d8a43-4233-4c83-8c1c-1c609ac23572
let
	l,m = -1,0.9375

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=60)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m09375.jpg", fig)
	
	fig
end

# ╔═╡ 4ef6016f-b0da-41cd-b880-f21bccb50a49
let
	l,m = -1,1

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 150
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=60)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m1.jpg", fig)
	
	fig
end

# ╔═╡ b8c6da2d-23e9-4db9-8942-edc7d6b1a513
let
	l,m = -1,1.9

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 80
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=10)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1_m19.jpg", fig)
	
	fig
end

# ╔═╡ 44aeb445-bd3f-4c0d-ba09-3f2900da684d
let
	l,m = -3ℯ/4,0.25

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=300)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-2038_m025.jpg", fig)
	
	fig
end

# ╔═╡ f7844fd7-3bd4-4cd7-886f-9a68c5eb4dd9
let
	l,m = -ℯ,0.0625

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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=100)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-2718_m00625.jpg", fig)
	
	fig
end

# ╔═╡ 27c14328-82b4-46e6-9906-b85a667e3dff
let
	l,m = -2.0794454,0.25

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=300)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-2079_m025.jpg", fig)
	
	fig
end

# ╔═╡ 7003b714-d25b-422a-9662-2bc56eca6aa8
let
	l,m = -ℯ,0.25

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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=100)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-2718_m025.jpg", fig)
	
	fig
end

# ╔═╡ 8e93509d-60b9-4f3a-8c30-8992fd42c73d
let
	l,m = -3.3,0.25

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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=0.5, markersize=4, iterations=200)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-33_m025.jpg", fig)
	
	fig
end

# ╔═╡ 3ac65d1b-4be0-4a15-9679-1dc485fde018
let
	l,m = -5ℯ/8,5

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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=2.5, markersize=8, iterations=3)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1))#=, Point2f(real(c1), -imag(c1))=#], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2))#=, Point2f(real(c2), -imag(c2))=#], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3))#=, Point2f(real(c3), -imag(c3))=#], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1698_m5.jpg", fig)
	
	fig
end

# ╔═╡ 1a785d2a-b284-465c-9585-a4a17388d608
let
	l,m = -9ℯ/8,4.75

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
		hasescaped = z -> real(z) < -32, 
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=2.5, markersize=8, iterations=5)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), #=Point2f(real(c1), -imag(c1))=#], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), #=Point2f(real(c2), -imag(c2))=#], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), #=Point2f(real(c3), -imag(c3))=#], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-3058_m475.jpg", fig)
	
	fig
end

# ╔═╡ 449ff397-4dbc-4163-b43f-2371e6c0b468
let
	l,m = -4ℯ/3,5.5

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
		hasescaped = z -> real(z) < -32, 
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=2.5, markersize=8, iterations=4)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), #=Point2f(real(c1), -imag(c1))=#], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), #=Point2f(real(c2), -imag(c2))=#], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), #=Point2f(real(c3), -imag(c3))=#], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-3624_m55.jpg", fig)
	
	fig
end

# ╔═╡ f51e3315-b30f-44d0-9ef3-23ec02c7bbdd
let
	l,m = -3ℯ/16,5.5

	Npix = 400
	xmin,xmax,ymin,ymax = -10,4,-12,12
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 100
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
		hasescaped = z -> real(z) < -12, 
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

	orbitpath!(ax, f, c1, colormap=[RGBA(1,0,0,0.5)], linewidth=1.5, markersize=4, iterations=200)
	#orbitpath!(ax, f2, conj(c1), colormap=[:yellow], linewidth=0.5, markersize=4, iterations=20)
	#orbitpath!(ax, f2, c2, colormap=[:magenta], linewidth=0.5, markersize=4, iterations=20)

	scatter!(ax, [Point2f(real(c1), imag(c1)), Point2f(real(c1), -imag(c1))], color=:red, markersize=8)	
	scatter!(ax, [Point2f(real(c2), imag(c2)), Point2f(real(c2), -imag(c2))], color=:blue, markersize=8)	
	scatter!(ax, [Point2f(real(c3), imag(c3)), Point2f(real(c3), -imag(c3))], color=:green, markersize=8)	

	text!(ax, Point2f(-7,0), text=L"F_\infty", color=:black, align=(:center,:center),
		  fontsize=0.75, markerspace=:data)

	Colorbar(fig[1, 2], limits = (0, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("DynPlane_l-1698_m5.jpg", fig)
	
	fig
end

# ╔═╡ 6b8acca5-4fd3-4812-a6c9-7f75f6b17b64
mandelbrot(
	(c,z)->z^2+c, -2.01:0.01:0.51,-1.16:0.01:1.16, seed=0, colormap=cuherx, 
	hasescaped = (c,z) -> abs2(z)>4, maxiterations=48,
	depth_shift=1, interpolate=true, fxaa=true, ssao=true, # transparency=false, alpha=0.9, 
	axis=(;aspect=DataAspect(), limits=(-2.01,0.51,-1.16,1.16), xticks=-2:0.25:0.5, yticks=-1:0.25:1, xgridcolor=RGBA(0.25,0.25,0.25,0.25), ygridcolor=RGBA(0.25,0.25,0.25,0.25), backgroundcolor=RGBA(1,1,1,0)))

# ╔═╡ 35911590-c861-465d-921c-de93b66de7e7
let
	Npix = 400
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
		hasescaped = (lm,z) -> abs2(z-f2(real(lm), imag(lm), z))<0.00001 && real(f(real(lm), imag(lm), z))<-16,
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

# ╔═╡ 8ac6cf60-315d-40ac-87fe-9af55c77f963
let
	Npix = 200
	xmin,xmax,ymin,ymax = -24,0,0,32
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuherx #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
	
	fig = Figure(size=(Npix,Npix))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:2:xmax, xgridcolor=gridc,
			yticks=ymin:2:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> abs2(z-f2(real(lm), imag(lm), z))<0.00001 && real(f(real(lm), imag(lm), z))<-16,
				#stops(z->f(real(lm), imag(lm), z), z, 2, ε=0.0000025) && real(f(real(lm), imag(lm), z))<-16, 
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	#lines!(ax,[-ℯ,0],[0,1],color=:red,linewidth=1.5)
	scatter!(ax, [Point2f(-20,0.25)], markersize=12, color=:red)
	text!(Point2f(-20,0.25), text="(-20,0.25)", color=:red, align = (:left,:bottom), fontsize=0.5,
		markerspace= :data)	

	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    label = "Iterations", vertical = true, flipaxis = true)

	#save("M1_plotbig.jpg", fig)
	
	fig
end

# ╔═╡ 0fa4f0d0-a16a-43d9-96ba-ab35321b0fe4
let
	Npix = 400
	xmin,xmax,ymin,ymax = -10,0,0,10
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 80
	cm = cuherx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,3Npix/4))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:1:xmax, xgridcolor=gridc,
			yticks=ymin:1:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticN(2, real(lm), imag(lm), maxiterations=100, maxiterationscompass=24, ε=0.000001),
		hasescaped = (lm,z) -> abs(real(z))>8 || abs2(z-f2(real(lm),imag(lm),z))<0.0000001,
		maxiterations=maxits, colormap=cm,
		interpolate=true, fxaa=true, ssao=true, depth_shift=1
	)

	fs = 0.096
	ms = 12
	
	#lines!(ax,[-ℯ,0],[0,1],color=:cadetblue,linewidth=1.75)
	#text!(Point2f(-0.18,0.74), text=L"\mathcal{L}_{-1}",
	#	  color=:dodgerblue4, align = (:center,:bottom), fontsize=fs, markerspace= :data, offset=(0.025,0.025))	
	
	Colorbar(fig[1, 2], limits = (1, maxits), colormap = cm,
    	label = "Iterations", vertical = true, flipaxis = true)	
	
	#save("E2_plot.jpg", fig)
	
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
	Npix = 200
	xmin,xmax,ymin,ymax = -3,0,0,1.5
	#xmin,xmax,ymin,ymax = -24,0,0,32
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 200
	cm = cuhex #vermeerx
	gridc = RGBA(0.25,0.25,0.25,0.25)
		
	fig = Figure(size=(Npix,Npix/2))
	ax = Makie.Axis(fig[1,1], aspect=DataAspect(), limits=(xmin,xmax,ymin,ymax),
			backgroundcolor=RGBA(1,1,1,0),
			xticks=xmin:2:xmax, xgridcolor=gridc,
			yticks=ymin:2:ymax, ygridcolor=gridc)

	mandelbrot!(ax,
		(lm,z) -> f2(real(lm), imag(lm), z), xs, ys,
		seed = lm -> findcriticBif(0.001,2.999, real(lm), imag(lm), maxiterations=32, ε=0.00001),
		hasescaped = (lm,z) -> stops(z->f(real(lm), imag(lm), z), z, 1, ε=0.0000025), 
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
	
	#save("BulbPer1Finite_plotbig.jpg", fig)
	
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
			backgroundcolor=RGBA(1,1,1,0.5),
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

# ╔═╡ 05a11fcd-afe9-4f9e-9e3d-b8d698085bf6
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 60
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
	
	#save("BulbBakerInf1_plot.jpg", fig)
	
	fig
end

# ╔═╡ cec5a4bd-a6ba-4fe9-be63-ed6322e0f8c1
let
	Npix = 400
	xmin,xmax,ymin,ymax = -6,0,0,6
	Δx,Δy = (xmax-xmin)/Npix, (ymax-ymin)/Npix
	Δ = min(Δx,Δy)
	xs = xmin:Δ:xmax
	ys = ymin:Δ:ymax

	maxits = 60
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
		hasescaped = (lm,z) -> real(z)<-16,
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
	
	#save("BulbBakerInf2_plot.jpg", fig)
	
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
	
	#save("Bifurcations_plot.pdf", fig)
	
	fig
end

# ╔═╡ 7e516588-5ad2-498b-89a6-91a376600ab2
cuherx_dark = pushfirst!(Gr.reverse(:cubehelix),RGB(0,0,0))

# ╔═╡ Cell order:
# ╟─afe66166-e9ca-11f0-bad5-59bd6959ad58
# ╠═dcc54357-7342-4978-a0e3-ae25b3161158
# ╠═04d2ed2f-de6c-4e5a-a5de-c848b6de8a87
# ╠═729aaaf9-28f3-4274-9a26-12fc6f80f074
# ╠═4662f2ab-8b41-4211-98e0-58778d61d007
# ╠═f09e37b5-b0af-46c3-bf25-51017bf48ffa
# ╠═a2cf3c8b-c5c8-4df5-951e-62a7954e1083
# ╠═9bceda1a-3a79-4e35-95b2-2572d5577993
# ╠═cada7f65-ccb1-4ec3-98e6-2370aac392cc
# ╟─72035530-a6b7-4c0e-8420-a4b2617d739e
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
# ╟─fdbde74c-6ede-422c-955f-de5e02ec8b22
# ╠═978cc14e-3383-450b-a46f-0376baa84b8c
# ╠═781de362-0f5f-44e5-87cc-d77345617ca7
# ╠═d2a17754-e606-48e9-82a1-a350091538d2
# ╠═0bce24f6-9278-464e-9d0c-0af471089dec
# ╠═e3130dbc-4b0d-4fa0-80a2-2980fd875d0b
# ╠═e8e876c1-a196-406a-b444-a8618a208690
# ╠═6b7870f5-a79d-416c-9a2d-6bd8a38e5fc5
# ╠═4cc91117-7c7b-4e70-820c-aff85860ab02
# ╟─208f363e-dc10-4307-aa27-3299ba2057a5
# ╟─e410e775-a417-40a2-8dbf-4c9c5050eef7
# ╟─70f58eb7-adc8-451c-94a6-11344bbb200b
# ╟─7273ed9e-5dfb-475b-b2a1-a65bfb64ab87
# ╟─aba0a028-b304-4a8d-91d3-f52693f34559
# ╟─12e25b26-5f8b-4bd8-9269-9b4672c66bf8
# ╠═0e4944e6-0626-4644-bcc2-6c3f1a5c5b37
# ╠═a8ede40d-b3b4-47f7-9a3c-2a787dc771a6
# ╟─028a0597-e08f-4424-b379-b7f604a0b79e
# ╟─9a344400-1172-4666-b4cd-1daec8f7339f
# ╠═3e83fb22-774d-4c41-85c6-21567d542ecc
# ╟─406d4a9f-539a-4510-aa56-6cb1e0dcfcf1
# ╟─01e1f581-df35-4179-9af2-3ce2e7a8ee20
# ╠═d188982b-663f-4405-916d-461ae5e24ee5
# ╟─fcb2bc42-4789-491b-b966-455569671de2
# ╟─5452962c-8c9b-4e27-bdb5-608259f16f49
# ╠═a983327e-8782-45d8-8459-a10dd7ad209b
# ╟─876d8a43-4233-4c83-8c1c-1c609ac23572
# ╟─e28fc392-446f-443e-97ed-fd7de84efd5f
# ╟─4ef6016f-b0da-41cd-b880-f21bccb50a49
# ╟─d1e50b30-e677-4fe8-8cb7-e65aac67441b
# ╟─263e7c61-b903-4cf7-9c8b-2872fb0c9989
# ╟─b8c6da2d-23e9-4db9-8942-edc7d6b1a513
# ╟─9e603559-2a3d-4d05-829a-99ef2dac32fb
# ╠═276212fc-720a-441d-ad15-0a3cc7c8c5a4
# ╟─44aeb445-bd3f-4c0d-ba09-3f2900da684d
# ╟─60e1123d-f204-4eca-a522-ebbc7d38e263
# ╠═a4fe1407-067c-4192-a6eb-844423a7a279
# ╟─f7844fd7-3bd4-4cd7-886f-9a68c5eb4dd9
# ╟─7bb26926-347f-4a43-8fa7-f3cff2bddfae
# ╟─27c14328-82b4-46e6-9906-b85a667e3dff
# ╟─d6c23140-9dba-49c8-8a64-4186ccf1ca9e
# ╠═795898a7-fabb-418f-aa19-ef98afb490f9
# ╟─7003b714-d25b-422a-9662-2bc56eca6aa8
# ╟─9bf7b086-d2d2-430a-a891-e43932c5c79d
# ╟─8e93509d-60b9-4f3a-8c30-8992fd42c73d
# ╟─7e8d748e-50ee-48a5-9081-d7f0e400647d
# ╠═9693d1da-7886-4d7a-906d-004196bd338f
# ╠═3ac65d1b-4be0-4a15-9679-1dc485fde018
# ╟─2bb4db46-4d9c-4eb0-8ea9-99d95d2804fa
# ╠═ac30ebe5-915b-4f2a-9634-912d4d7ad369
# ╟─1a785d2a-b284-465c-9585-a4a17388d608
# ╟─e17ec7eb-b95e-4e93-894c-3dab3cf26e10
# ╠═9f02bbfc-d1c4-43d7-a940-a6b6a6b7305b
# ╟─449ff397-4dbc-4163-b43f-2371e6c0b468
# ╟─f51e3315-b30f-44d0-9ef3-23ec02c7bbdd
# ╟─0dc154e3-b602-42be-9b56-3325f1961581
# ╠═369f53e9-4ea9-4462-8297-ebb3fc381ffb
# ╠═edf8fc8f-b6bd-4f29-af08-00039b30acd1
# ╠═70a90845-4d5f-44e0-997b-4ff969b56346
# ╟─ffc2b1aa-c1d9-4c69-bccd-94d24b9744aa
# ╟─cdf5ac1b-df7f-4270-8c67-1340e20861b4
# ╠═99dc3a03-9aa0-4f52-8ed7-ffb4d4e5b538
# ╠═7786dcfb-1d6e-419b-802c-e812a37361eb
# ╟─0fd2c626-52e9-4f8c-8eac-64a2047d9729
# ╟─6b8acca5-4fd3-4812-a6c9-7f75f6b17b64
# ╟─d91916f8-37a6-4494-ac61-6df7aba6a0c1
# ╟─6e0bf940-1d94-4efa-a0a0-6558ba5d67fe
# ╠═4c252064-eb11-485b-9607-01aaf84242bb
# ╟─735528fa-e490-47ab-bc6c-bbba9b2fdac8
# ╠═0ec76c85-d021-4a0b-8872-882c03f1a0d3
# ╠═35911590-c861-465d-921c-de93b66de7e7
# ╠═8ac6cf60-315d-40ac-87fe-9af55c77f963
# ╟─65365633-6555-4f3d-945b-66275dcf8c7b
# ╟─fef962ec-4550-4a2f-bb28-5df9c7ced028
# ╟─ffbacbc1-1b56-4b1f-93f1-0dce2f44bc9f
# ╠═74aa3553-400b-4df5-b053-7ee27338d8b6
# ╟─249dda01-d980-4cc6-86f0-1df2a31dd7d6
# ╟─07765e47-4d63-4df8-ab8f-906d71a75338
# ╟─2b61308c-3efa-47a9-88d8-0ea870ca3d77
# ╟─d2598cd0-8ede-4558-bc38-b5503bab5f17
# ╠═2176a612-f784-441b-8d44-ba132d5883b8
# ╠═855b2724-3129-4166-b170-99f6ae6cbc37
# ╠═a46480b9-1b34-4c81-82e1-be4af1ed3e05
# ╟─ab526193-7dfd-4e82-a397-beb7e2fd6c93
# ╟─e23f79d7-5899-409b-b721-4bb06aac37b4
# ╟─08e7bb0a-e954-46e0-9be8-eec3fce7d864
# ╠═05a11fcd-afe9-4f9e-9e3d-b8d698085bf6
# ╠═cec5a4bd-a6ba-4fe9-be63-ed6322e0f8c1
# ╟─40f1cf47-2b5c-4d9e-bc8b-10ca781dffcc
# ╠═3d076c9b-263b-4d62-a17e-afbd7fee2f21
# ╟─05efe4d3-4126-449a-ac48-1799c4c2f28a
# ╟─de1c003f-d951-47cf-acc3-beaea021298b
# ╟─2c9f1a59-ad76-4476-8d32-7e158dcd2541
# ╠═ee95a23b-177b-474d-a517-5cc9204831dd
# ╠═0fa4f0d0-a16a-43d9-96ba-ab35321b0fe4
# ╟─5cc6c94c-fccf-4251-aaf3-0629f396457f
# ╠═fab68d9a-311d-451f-9289-caaec83569dc
# ╠═d8613d71-af5e-4eb9-b065-b6467eae28b0
# ╟─ee2b8e91-d405-4045-9287-e72daf20c7fd
# ╠═cd230555-d920-4eb7-863c-e4ef586e84a0
# ╠═d52d0e77-9a16-4163-b2ef-10779cd1cdff
# ╠═7a0d58fd-7238-4f0c-accf-9e429ca4d325
# ╠═7af6c311-4ff3-4f71-84f0-8d9ef2d591c8
# ╠═df10e677-8b09-473e-8b15-3ee063496f34
# ╠═92366d0e-43c7-49c9-b78a-07ced091dcee
# ╠═36239fa1-7865-4b54-bd6a-025c80bceae8
# ╠═4b8a3a21-889f-4189-a6c5-e06014e1e8e0
# ╠═4294274f-8874-4840-a241-0cce4969d5dd
# ╠═a07e4983-d801-47ba-9531-16bdfb1f5f26
# ╠═c14eacb3-f173-4609-b750-413027a296ca
# ╠═fd700d57-b159-4a0e-87e6-6b7361254b60
# ╠═57f7c4ba-d10a-4f0d-ae0d-b12f49a1c313
# ╠═0947a4b6-0309-4fea-bcab-7a3f0d036974
# ╠═7e516588-5ad2-498b-89a6-91a376600ab2
