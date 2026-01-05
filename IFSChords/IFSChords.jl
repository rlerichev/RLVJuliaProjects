### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ f72af3dc-54d7-40bf-97c2-ae0621c77ef0
begin
	using Pkg
	urlsdd="https://github.com/Colectivo-SDD/SDD"
	Pkg.add(url=urlsdd*"Core.jl")
	Pkg.add(url=urlsdd*"Graphics.jl")
	Pkg.add(url=urlsdd*"Geometry.jl")
	Pkg.add(url=urlsdd*".jl") #, rev="rlv_makie")	
	#Pkg.add(url=urlsdd*"IFS.jl", rev="rlv_makie")
end

# ╔═╡ 463d3d7a-334e-44e1-9369-9422ad8c48aa
using Colors, ColorSchemes, Images

# ╔═╡ 260ae4ae-c2cf-42dd-88f3-860435d1da88
using SDDCore, SDDGeometry, SDDGraphics, SDD

# ╔═╡ 50981a90-f37c-47ad-9d2b-abd1bd0af77e
md"""
## Función continua de las inversas con pesos
Se construirá una función continua, $F:\mathbb{C}\rightarrow\mathbb{C}$, cuyo conjunto de puntos atrapados representa información musical en un tiempo dado (*frame*).

---

Sean $\{p_1,\dots,p_N\}$ los puntos en el (los) círculo(s) cromático(s) que representan el tono (*pitch*).

Las funciones contractivas afines para el atractor fractal son de la forma

$f_k(z)=a_k(z-p_k)+p_k,$

donde $a_k\in\mathbb{D}-\{0\}$.

Entonces las correspondientes inversas son

$g_k(z)=\frac{1}{a_k}(z-p_k)+p_k.$

---

Queremos que $F(p_k)=g_k(p_k)=p_k$, y para puntos $z$ cercanos a $p_k$ se tenga $F(z)\approx g_k(z)$. Una función así puede ser

$F(z)=\sum_{k=1}^N w_k(z)g_k(z)$

donde $w_k:\mathbb{C}\rightarrow\mathbb{C}$ son funciones peso (*weight*) tales que $w_k(p_k)=1$ y $w_k(p_j)=0$ si $k\neq j$.

---

Funciones como las $w_k$ son los **polinomios de Lagrange** $\ell_k$. Dados $\{p_1,\dots,p_N\}$ y un $k\in\{1,\dots,N\}$, se define

$\ell_k(z)=\prod_{j=1,j\neq k}^N \frac{z-p_j}{p_k-p_j}=\frac{(z-p_1)\dots(z-p_{k-1})(z-p_{k+1})\dots(z-p_N)}{(p_k-p_1)\dots(p_k-p_{k-1})(p_k-p_{k+1})\dots(p_k-p_N)}.$


Forma baricéntrica, para cómputo eficiente:

$\ell(z)=\prod_{j=1}^N (z-p_j),$

$b_k=\prod_{j=1,j\neq k}^N \frac{1}{p_k-p_j},$

$\ell_k(z)=\ell(z)\frac{b_k}{z-p_k}=
\begin{cases}
1 && z=p_k, \\
0 && z=p_j,\,j\neq k,\\
\frac{b_k}{(z-p_k)\sum_{m=1}^N\frac{b_m}{z-p_m}} && z\neq p_j.
\end{cases}$

La derivada:

$\ell'_k(z)=\ell_k(z)\sum_{j=1,j\neq k}^N \frac{1}{z-p_j}.$

$\implies$

$\ell'_k(p_k)=\ell_k(p_k)\sum_{j=1,j\neq k}^N \frac{1}{p_k-p_j}=\sum_{j=1,j\neq k}^N \frac{1}{p_k-p_j},$

$\ell'_k(p_m)=\ell_k(p_m)\sum_{j=1,j\neq k}^N \frac{1}{p_m-p_j}=\sum_{j=1,j\neq k,m}^N \frac{p_k-p_m}{p_k-p_j}???.$



---
Nota: Checar el paquete **SpecialPolynomials.jl**, que ya implementa polinomios de Lagrange...

"""

# ╔═╡ 5d3a715c-619a-42b6-bf47-11e84744c4e7
md"""
---

### 1

Primera propuesta de funciones peso:

$w_k(z)=e^{-|z-p_k|^2}|\ell_k(z)|^2.$

Nótese que
-  $F$ sería diferenciable pero no holomorfa.


-  $w_k:\mathbb{C}\rightarrow\mathbb{R}$.


Otras posibilidades parecidas a esta propuesta son $w_k(z)=|\ell_k(z)|$, $w_k(z)=|\ell_k(z)|^2$, $w_k(z)=|\frac{z}{p_k}\ell_k(z)|$, etc.

"""

# ╔═╡ cfe72a01-6b0e-45f7-84fe-36a421906730
md"""
---

### 2

Segunda propuesta de funciones peso:

$w_k(z)=\frac{z}{p_k}\ell_k(z).$

Nótese que
-  $F$ sería un polinomio y por lo tanto holomorfa. Al ser polinomio se tiene que $\infty$ es super-atractor y entonces el conjunto de puntos atrapados es un compacto en $\mathbb{C}$, por lo que este conjunto siempre se puede dibujar bien con el algoritmo de coloreado por tiempo de escape.


-  $w_k:\mathbb{C}\rightarrow\mathbb{C}$.


-  $F(0)=0$.

"""

# ╔═╡ 7634247e-b2ac-489f-943a-cef75f5a512b
md"""
---

### 3

Tercera propuesta de funciones peso:

$w_k(z)=\ell_k(z).$


Nótese que
-  $F$ sería un polinomio y por lo tanto holomorfa.


-  $w_k:\mathbb{C}\rightarrow\mathbb{C}$.


-  $F(0)\neq 0$ en general.

"""

# ╔═╡ b9a44483-41c8-41df-8ade-3283c23676e4
md"""
---

### 4

Cuarta propuesta de funciones peso (ESTÁN MAL ESTAS FUNCIONES, NO SE CUMPLE LO ESPERADO DE LA DERIVADA):

$w_k(z)=\ell_k(z)-(z-p_k)\ell'_k(p_k).$


Nótese que
-  $F$ sería un polinomio y por lo tanto holomorfa.


-  $w_k:\mathbb{C}\rightarrow\mathbb{C}$.


-  $F(0)\neq 0$ en general.


-  $F'(p_k)=\frac{1}{a_k}$, que nos asegura que $p_k$ es fijo repulsor (es lo que quisiéramos).
"""

# ╔═╡ 8572c6aa-524d-4024-b6e5-2ec7a77b51d5
md"""
---

Algunos cálculos, tal vez útiles...

$F'(z)=\sum_{k=1}^N (w'_k(z)g_k(z)+w_k(z)g'_k(z))=\sum_{k=1}^N(w'_k(z)g_k(z)+w_k(z)\frac{1}{a_k})$

$F'(p_m)= \sum_{k=1}^N (w'_k(p_m)g_k(p_m)+w_k(p_m)\frac{1}{a_k})$

$=w'_m(p_m)p_m+\frac{1}{a_m}+\sum_{k=1,k\neq m}^N w'_k(p_m)g_k(p_m)$

"""

# ╔═╡ cfd5ba6c-0af9-4613-a5ec-d97d2320f749
md"""
---

### 5

Cuencas de atracción.

Tomando la función del método de Newton-Raphson de la función cuyas raíces son los puntos fijos de $F$, es decir, los puntos asociados a los tonos $\{p_1,\dots,p_N\}$.

$G(z)=F(z)-z.$

$N_G(z)=z-\frac{G(z)}{G'(z)}.$


"""

# ╔═╡ 86aa468a-5091-4f24-a600-b46485709d44
md"""
### 6

Exponencial con racional.

$w_k(z)=e^{-\Big(\frac{z-p_k}{(z-p_1)\dots(z-p_{k-1})(z-p_{k+1})\dots(z-p_N)}\Big)}.$

(No está completamente bien, hay que asegurar que para $j\neq k$ se obtenga $e^{-\infty_{\mathbb{R}}}=0$.)
"""

# ╔═╡ dd76ba61-2570-47ad-9a04-21af3f6b1ade
R(z::Number) = (z-2-im)/(z-1+im)

# ╔═╡ 59d391b6-28ba-4ecc-9750-b7c44c2d327e
eR(z) = exp(-R(z))

# ╔═╡ 176ae08d-7f41-4d31-8a03-eaca808bd682
eR(1-im)

# ╔═╡ 65e63390-b2c4-4b85-aec6-7e0f81bcff3a
md"""
Exponencial con racional no holomorfa.

$w_k(z)=e^{-\Big|\frac{z-p_k}{(z-p_1)\dots(z-p_{k-1})(z-p_{k+1})\dots(z-p_N)}\Big|}.$

Estas si cumplen que $w_k(p_j)=0$.
"""

# ╔═╡ 1b4ebdf8-e7b5-427a-a69e-3cbd5a6de6f3
md"""
# Ejemplos

En los siguientes ejemplos los colores no tienen interpretación musical. Sólo se muestran las "formas esenciales" de algunos acordes.
"""

# ╔═╡ 767e88c6-da17-11f0-abd7-dfe2ab96d6b7
md"""
## Tricordio
"""

# ╔═╡ 735dc78a-f829-4a55-b7bb-abbbf15813fb
function f3chord(z::Number, p1::Number, p2::Number, p3::Number, s::Number=8)
	r1 = exp(-abs(z-p1))*abs((z-p2)*(z-p3)/((p1-p2)*(p1-p3)))
	r2 = exp(-abs(z-p2))*abs((z-p3)*(z-p1)/((p2-p3)*(p2-p1)))
	r3 = exp(-abs(z-p3))*abs((z-p1)*(z-p2)/((p3-p1)*(p3-p2)))
	#r1 = exp(-abs((z-p1)/((z-p2)*(z-p3))))
	#r2 = exp(-abs((z-p2)/((z-p3)*(z-p1))))
	#r3 = exp(-abs((z-p3)/((z-p1)*(z-p2))))
	r1*(s*(z-p1)+p1) + r2*(s*(z-p2)+p2) + r3*(s*(z-p3)+p3)
end

# ╔═╡ 95ecba59-cbba-4f09-b819-2cdc476f39df
md"""
---

Do mayor (C), primera propuesta
"""

# ╔═╡ ade51186-aa55-4250-bddb-f859e634c987
imgtrappedpoints(z->f3chord(z,1,exp(2pi*im*4/12),exp(2pi*im*7/12),4.25), -6:0.01:6, -6:0.01:6, colormap=:jet, maxiterations=20, hasescaped=z->abs2(z)>36)

# ╔═╡ d547d1fa-0148-40ce-bc6a-503fdda95f3f
md"""
---

C tritono (C+5?), primera propuesta
"""

# ╔═╡ ed704397-87af-43a8-b1ca-989295504768
imgtrappedpoints(z->f3chord(z,1,exp(2pi*im*4/12),exp(2pi*im*8/12),4.25), -6:0.01:6, -6:0.01:6, colormap=:jet, maxiterations=20, hasescaped=z->abs2(z)>36)

# ╔═╡ dc957f23-7083-40ae-8e34-337c748912f5
md"""
---
"""

# ╔═╡ 4eb5a094-a3d4-4725-9dfe-18b0bd506c68
function f3chordpoly(z::Number, p1::Number, p2::Number, p3::Number, s::Number=2.5)
	r1 = z*(z-p2)*(z-p3)/(p1*(p1-p2)*(p1-p3))
	r2 = z*(z-p3)*(z-p1)/(p2*(p2-p3)*(p2-p1))
	r3 = z*(z-p1)*(z-p2)/(p3*(p3-p1)*(p3-p2))
	r1*(s*(z-p1)+p1) + r2*(s*(z-p2)+p2) + r3*(s*(z-p3)+p3)
end

# ╔═╡ 81830da4-fcb4-488d-bf95-c6444a9b1ed4
md"""
---

Do mayor (C), segunda propuesta
"""

# ╔═╡ 694bd4db-4f78-4ee1-89b5-9001cf2a6ee1
imgtrappedpoints(z->f3chordpoly(z,1,exp(2pi*im*3/12),exp(2pi*im*7/12)), -1.25:0.005:1.25, -1.25:0.005:1.25, colormap=:jet, maxiterations=24, hasescaped=z->abs2(z)>4)

# ╔═╡ c2f3822d-dbec-48e7-bfb2-d617693fa7f5
f3chordpoly(0,1,exp(2pi*im*3/12),exp(2pi*im*7/12))

# ╔═╡ 1ff50618-ca25-40f0-9e1f-0c787713c9c1
md"""
---

C tritono (C+5?), primera propuesta
"""

# ╔═╡ 50366e11-44cf-4158-aa0b-008b9c69836b
imgtrappedpoints(z->f3chordpoly(z,1,exp(2pi*im*4/12),exp(2pi*im*8/12),2.3), -1.5:0.005:1.5, -1.5:0.005:1.5, colormap=:jet, maxiterations=20, hasescaped=z->abs2(z)>4)

# ╔═╡ 0c7a0502-5662-4e32-a69a-e9d14dc14660
md"""
---
"""

# ╔═╡ fabbcc12-5a00-4825-96b8-471757c81969
md"""
Tercera propuesta.

En realidad, se está usando $(w_k(z))^2$, pues sin el cuadrado no funciona bien...
"""

# ╔═╡ ef7936e8-d774-46ce-aac7-ddd34d5d4a77
function f3chordpoly3(z::Number, p1::Number, p2::Number, p3::Number, s::Number=2.5)
	r1 = ((z-p2)*(z-p3)/((p1-p2)*(p1-p3)))^2
	r2 = ((z-p3)*(z-p1)/((p2-p3)*(p2-p1)))^2
	r3 = ((z-p1)*(z-p2)/((p3-p1)*(p3-p2)))^2
	r1*(s*(z-p1)+p1) + r2*(s*(z-p2)+p2) + r3*(s*(z-p3)+p3)
end

# ╔═╡ 2e4d349c-1dde-4e08-bc90-3f3bb5f7cf64
imgtrappedpoints(z->f3chordpoly3(z,1,exp(2pi*im*4/12),exp(2pi*im*7/12),5.75), -1.5:0.005:1.5, -1.5:0.005:1.5, colormap=:jet, maxiterations=40, hasescaped=z->abs2(z)>24)

# ╔═╡ f875a301-8f2a-4dc8-aaa5-14f97176d74a
f3chordpoly3(exp(2pi*im*7/12),1,exp(2pi*im*4/12),exp(2pi*im*7/12),1.1)

# ╔═╡ 901619a7-b4af-44bf-8bb3-3a824b83a687
function f3chordpoly31(z::Number, p1::Number, p2::Number, p3::Number, s::Number=2.5)
	r1 = ((z-p2)*(z-p3)/((p1-p2)*(p1-p3)))
	r2 = ((z-p3)*(z-p1)/((p2-p3)*(p2-p1)))
	r3 = ((z-p1)*(z-p2)/((p3-p1)*(p3-p2)))
	r1*(s*(z-p1)+p1) + r2*(s*(z-p2)+p2) + r3*(s*(z-p3)+p3)
end

# ╔═╡ 9a38ff27-3fb2-4729-b1fb-2dab0a3b4dba
md"""
---

Cuencas de atracción
"""

# ╔═╡ b43e6a23-36f2-4e54-8fbf-a8ce1c426679
function numder(f::Function, δ::Real=0.0001)
  function fder(z::Number)
    (f(z+δ)-f(z-δ) - im*(f(z+δ*im)-f(z-δ*im)))/(4δ)
  end
end

# ╔═╡ a88e7c29-d196-4688-a431-63020fbf98a6
numder(sin)(0)

# ╔═╡ 1d730b36-0018-4ea8-b960-572bb6453782
function NewtonRaphsonFunc(f::Function, δ::Real=0.0001)
  function NRF(z::Number)
	  f´z = numder(f,δ)(z)
	  if abs2(f´z) < δ
		  return Inf
	  end
	  z - f(z)/f´z
  end
end

# ╔═╡ 0d2a1420-5529-47ef-ad72-e94c25fdf75b
imgbasins(
	NewtonRaphsonFunc(
		z->f3chordpoly(z,1,exp(2pi*im*4/12),exp(2pi*im*7/12),2)-z
	),
	[0,1,exp(2pi*im*4/12),exp(2pi*im*7/12)], -2:0.0075:2, -2:0.0075:2, colormap=:jet, maxiterations=16, tolerance=0.01)

# ╔═╡ d46afa0c-61ea-4ef1-a92c-b02c252514b0
imgbasins(
	NewtonRaphsonFunc(
		z->z^3-1, 0.000001
	),
	[1,exp(2pi*im/3),exp(4pi*im/3)], -2:0.01:2, -2:0.01:2, colormap=[RGB(1,1,1),RGB(0.75,0,0), RGB(1,1,0), RGB(0,0.25,0.5)], maxiterations = 24, tolerance=0.05)

# ╔═╡ 689f52af-a2e1-4c67-9ae8-d6fc435810db
imgbasins(
	z->f3chordpoly(z,1,exp(2pi*im*4/12),exp(2pi*im*7/12),0.5),
	[0,1,exp(2pi*im*4/12),exp(2pi*im*7/12)], -2:0.0075:2, -2:0.0075:2, colormap=:jet, maxiterations=16, tolerance=0.001)

# ╔═╡ a4549efd-ccf9-4d45-a1e6-70807e779f2e
md"""
## Tetracordio
"""

# ╔═╡ 2a83488a-c980-4784-bb23-5b881bb0f60d
function f4chord(z::Number, p1::Number, p2::Number, p3::Number, p4::Number, s::Real=10)
	r1 = exp(-abs2(z-p1))*abs2((z-p2)*(z-p3)*(z-p4)/((p1-p2)*(p1-p3)*(p1-p4)))
	r2 = exp(-abs2(z-p2))*abs2((z-p3)*(z-p4)*(z-p1)/((p2-p3)*(p2-p4)*(p2-p1)))
	r3 = exp(-abs2(z-p3))*abs2((z-p4)*(z-p1)*(z-p2)/((p3-p4)*(p3-p1)*(p3-p2)))
	r4 = exp(-abs2(z-p4))*abs2((z-p1)*(z-p2)*(z-p3)/((p4-p1)*(p4-p2)*(p4-p3)))
	r1*(s*(z-p1)+p1) + r2*(s*(z-p2)+p2) + r3*(s*(z-p3)+p3) + r4*(s*(z-p4)+p4)
end

# ╔═╡ 317a6d1f-3efd-4747-958a-1f01f5ec9a90
md"""
---

Do mayor con séptima (C7), primera propuesta.
"""

# ╔═╡ c69a7e92-e707-41b2-9119-9ed5366f5615
imgtrappedpoints(z->f4chord(z,1,exp(2pi*im*4/12),exp(2pi*im*7/12),exp(2pi*im*10/12),10.6), -6:0.01:6, -6:0.01:6, colormap=:jet, maxiterations=20, hasescaped=z->abs2(z)>36)

# ╔═╡ 2e49103c-61a8-4ab4-81e2-e098927d9856
md"""
---

Do disminuido (Cdim), primera propuesta.
"""

# ╔═╡ e2979324-dd3a-46ee-9de2-beca02fc6bb7
imgtrappedpoints(z->f4chord(z,1,exp(2pi*im*3/12),exp(2pi*im*6/12),exp(2pi*im*9/12),16), -6:0.01:6, -6:0.01:6, colormap=:jet, maxiterations=20, hasescaped=z->abs2(z)>36)

# ╔═╡ 94a441e1-a04c-489a-a884-337408703297
function f4chordpoly(z::Number, p1::Number, p2::Number, p3::Number, p4::Number, s::Real=2.5)
	r1 = z*(z-p2)*(z-p3)*(z-p4)/(p1*(p1-p2)*(p1-p3)*(p1-p4))
	r2 = z*(z-p3)*(z-p4)*(z-p1)/(p2*(p2-p3)*(p2-p4)*(p2-p1))
	r3 = z*(z-p4)*(z-p1)*(z-p2)/(p3*(p3-p4)*(p3-p1)*(p3-p2))
	r4 = z*(z-p1)*(z-p2)*(z-p3)/(p4*(p4-p1)*(p4-p2)*(p4-p3))
	r1*(s*(z-p1)+p1) + r2*(s*(z-p2)+p2) + r3*(s*(z-p3)+p3) + r4*(s*(z-p4)+p4)
end

# ╔═╡ 46b145eb-36f6-4ab2-98b0-d97ad7536108
md"""
---

Do mayor con séptima (C7), segunda propuesta.
"""

# ╔═╡ e62913a0-a415-4823-9f43-540a9d90636f
imgtrappedpoints(z->f4chordpoly(z, 1, exp(2pi*im*4/12), exp(2pi*im*7/12), exp(2pi*im*10/12),2.46), -1.5:0.005:1.5, -1.5:0.005:1.5, colormap=:jet, maxiterations=24, hasescaped=z->abs2(z)>4)

# ╔═╡ bb456fc2-e1b8-4dbf-a368-86d4c656b396
md"""
---

Do disminuido (Cdim), segunda propuesta.
"""

# ╔═╡ 9fc2814f-1e38-4289-aef5-b6bbb6081798
imgtrappedpoints(z->f4chordpoly(z, 1, exp(2pi*im*3/12), exp(2pi*im*6/12), exp(2pi*im*9/12),2.65), -1.5:0.005:1.5, -1.5:0.005:1.5, colormap=:jet, maxiterations=24, hasescaped=z->abs2(z)>4)

# ╔═╡ 7c022824-c3e8-41b5-be5a-458d56650cff
md"""
### Bibliotecas
"""

# ╔═╡ 1d27accf-c430-4435-9bab-3b0a4bb2d015
const Gr = SDDGraphics

# ╔═╡ Cell order:
# ╟─50981a90-f37c-47ad-9d2b-abd1bd0af77e
# ╟─5d3a715c-619a-42b6-bf47-11e84744c4e7
# ╟─cfe72a01-6b0e-45f7-84fe-36a421906730
# ╟─7634247e-b2ac-489f-943a-cef75f5a512b
# ╟─b9a44483-41c8-41df-8ade-3283c23676e4
# ╟─8572c6aa-524d-4024-b6e5-2ec7a77b51d5
# ╟─cfd5ba6c-0af9-4613-a5ec-d97d2320f749
# ╟─86aa468a-5091-4f24-a600-b46485709d44
# ╠═dd76ba61-2570-47ad-9a04-21af3f6b1ade
# ╠═59d391b6-28ba-4ecc-9750-b7c44c2d327e
# ╠═176ae08d-7f41-4d31-8a03-eaca808bd682
# ╟─65e63390-b2c4-4b85-aec6-7e0f81bcff3a
# ╟─1b4ebdf8-e7b5-427a-a69e-3cbd5a6de6f3
# ╟─767e88c6-da17-11f0-abd7-dfe2ab96d6b7
# ╠═735dc78a-f829-4a55-b7bb-abbbf15813fb
# ╟─95ecba59-cbba-4f09-b819-2cdc476f39df
# ╠═ade51186-aa55-4250-bddb-f859e634c987
# ╟─d547d1fa-0148-40ce-bc6a-503fdda95f3f
# ╠═ed704397-87af-43a8-b1ca-989295504768
# ╟─dc957f23-7083-40ae-8e34-337c748912f5
# ╠═4eb5a094-a3d4-4725-9dfe-18b0bd506c68
# ╟─81830da4-fcb4-488d-bf95-c6444a9b1ed4
# ╠═694bd4db-4f78-4ee1-89b5-9001cf2a6ee1
# ╠═c2f3822d-dbec-48e7-bfb2-d617693fa7f5
# ╟─1ff50618-ca25-40f0-9e1f-0c787713c9c1
# ╠═50366e11-44cf-4158-aa0b-008b9c69836b
# ╟─0c7a0502-5662-4e32-a69a-e9d14dc14660
# ╟─fabbcc12-5a00-4825-96b8-471757c81969
# ╠═ef7936e8-d774-46ce-aac7-ddd34d5d4a77
# ╠═2e4d349c-1dde-4e08-bc90-3f3bb5f7cf64
# ╠═f875a301-8f2a-4dc8-aaa5-14f97176d74a
# ╠═901619a7-b4af-44bf-8bb3-3a824b83a687
# ╠═9a38ff27-3fb2-4729-b1fb-2dab0a3b4dba
# ╠═b43e6a23-36f2-4e54-8fbf-a8ce1c426679
# ╠═a88e7c29-d196-4688-a431-63020fbf98a6
# ╠═1d730b36-0018-4ea8-b960-572bb6453782
# ╠═0d2a1420-5529-47ef-ad72-e94c25fdf75b
# ╠═d46afa0c-61ea-4ef1-a92c-b02c252514b0
# ╠═689f52af-a2e1-4c67-9ae8-d6fc435810db
# ╟─a4549efd-ccf9-4d45-a1e6-70807e779f2e
# ╠═2a83488a-c980-4784-bb23-5b881bb0f60d
# ╟─317a6d1f-3efd-4747-958a-1f01f5ec9a90
# ╠═c69a7e92-e707-41b2-9119-9ed5366f5615
# ╟─2e49103c-61a8-4ab4-81e2-e098927d9856
# ╠═e2979324-dd3a-46ee-9de2-beca02fc6bb7
# ╠═94a441e1-a04c-489a-a884-337408703297
# ╟─46b145eb-36f6-4ab2-98b0-d97ad7536108
# ╠═e62913a0-a415-4823-9f43-540a9d90636f
# ╟─bb456fc2-e1b8-4dbf-a368-86d4c656b396
# ╠═9fc2814f-1e38-4289-aef5-b6bbb6081798
# ╟─7c022824-c3e8-41b5-be5a-458d56650cff
# ╠═463d3d7a-334e-44e1-9369-9422ad8c48aa
# ╠═f72af3dc-54d7-40bf-97c2-ae0621c77ef0
# ╠═260ae4ae-c2cf-42dd-88f3-860435d1da88
# ╠═1d27accf-c430-4435-9bab-3b0a4bb2d015
