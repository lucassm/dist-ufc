### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3add8c37-e20b-48ca-b6b7-db6d3e38e6f9
begin
	using LinearAlgebra
	using Statistics
	using Printf
	using PlutoUI
end

# ╔═╡ deb8bd3d-cffd-49af-8d16-35122a923ebc
md"""
# Fluxo de Potência com Transformador Trifásico D-Y aterrado  
"""

# ╔═╡ 3dcaa6fe-080a-4fb1-b18c-871249e513d6
md"""
## Example 8.2

In the example system in Figure, an unbalanced load is being served at the end of a **10,000-ft** section of a three-phase line.

> The 10,000 ft long line is being fed from a substation transformer rated 5000 kVA, 115 kV delta—12.47 kV grounded wye with a per-unit impedance of 0.085 /_ 85º.

> The phase conductors of the line are 336,400 26/7 Aluminum Conductor Steel Reinforced (ACSR) with a neutral conductor 4/0 ACSR.

The configuration and computation of the phase impedance matrix are given in Example 4.1. From that example, the phase impedance matrix was computed to be:

$\left[z_{\text {line }}\right]=\left[\begin{array}{lll}
0.4576+j 1.0780 & 0.1560+j 0.5017 & 0.1535+j 0.3849 \\
0.1560+j 0.5017 & 0.4666+j 1.0482 & 0.1580+j 0.4236 \\
0.1535+j 0.3849 & 0.1580+j 0.4236 & 0.4615+j 1.0651
\end{array}\right] \Omega / \mathrm{mile}$

The source voltages at node 1 are:

$\left[E L L_{A B C}\right]=\left[\begin{array}{c}
115,000 \angle 0^\circ \\
115,000 \angle -120^\circ \\
115,000 \angle 120^\circ
\end{array}\right]$

The wye-connected loads are:

$[k V A]=\left[\begin{array}{l}
1700 \\
1200 \\
1500
\end{array}\right] \quad[P F]=\left[\begin{array}{l}
0.90 \\
0.85 \\
0.95
\end{array}\right]$

The complex powers of the loads are computed to be:

$SL_i
=
k V A_i \cdot \mathrm{e}^{j \cdot \operatorname{acos}(P F i)}
=
\left[\begin{array}{l}
1530+j 741.0 \\
1020+j 632.1 \\
1425+j 468.4
\end{array}\right]
\mathrm{kW}+j k v a r$

"""

# ╔═╡ d36ab22c-2897-4dad-a693-c70eaa1f81dc
PlutoUI.TableOfContents(include_definitions=false)

# ╔═╡ 2ff6e94d-872a-471a-9ce3-dfe2a491f31e
md"""
## Montagem da Matriz de Impedâncias Primitiva
"""

# ╔═╡ 16bcfca2-557b-4512-a4b5-7e05c98382cd
D = [0.0244 2.5 7.0 5.6569;
	 2.5 0.0244 4.5 4.272;
	 7.0 4.5 0.0244 5.0;
	 5.6569 4.272 5.0 0.0081]

# ╔═╡ 623e3772-f53a-4765-a365-1117aa5110ae
begin
	Zabc = zeros(Complex, 4, 4)
	for i in 1:4
		for j in 1:4
			if i != j
				Zabc[i, j] = 0.0953 + 0.12134 * (log(1.0 / D[i, j]) + 7.93402)im
			elseif i == 4
				Zabc[i, j] = 0.592 + 0.0953 + 0.12134 * (log(1.0 / D[i, j]) + 7.93402)im
			else
				Zabc[i, j] = 0.306 + 0.0953 + 0.12134 * (log(1.0 / D[i, j]) + 7.93402)im
			end
		end
	end
end

# ╔═╡ 960141d1-7b14-4e05-8e6a-aa37484d2b8b
Zabc

# ╔═╡ 664eb69b-9fc8-4c98-a899-40aac945d3de
md"""
## Montagem da Matriz de Impedâncias de Fase (Redução de Kron)
"""

# ╔═╡ d4c98edb-74ff-46ab-b8b0-d7acb115cfc7
md"""
Redução de Kron da matriz de impedâncias primitiva
"""

# ╔═╡ d6b1dc9c-c3a9-46ea-ac64-937f31f35e4e
Zabck = Zabc[1:3, 1:3] - Zabc[1:3, 4:4] * Zabc[4:4, 1:3] / Zabc[4, 4]

# ╔═╡ 248d7de5-6b58-495b-a788-2450706c5b9f
l = 10e3 / 5280.0

# ╔═╡ 478ddb02-083d-45a0-9aad-6459235a878a
md"""
Matriz de Impedâncias Primitiva em Ohms:
"""

# ╔═╡ 571affb9-0c8f-4541-a4b1-3bea3e45b036
Zabc_ = Zabc * l

# ╔═╡ 542462a7-8c53-422d-96ba-d4ce6db0deda
md"""
Matriz de Impedâncias de Fase Reduzida de Kron em Ohms:
"""

# ╔═╡ 2de5a128-85b8-4111-9dd5-5c6450d7ad1b
Zabck_ = Zabck * l

# ╔═╡ e45f3e9f-71c6-4a63-90a0-ee635dec4a43
md"""
## Definição das Matrizes características da Linha
"""

# ╔═╡ a51756e4-d14f-430f-b18b-5beda9f321cf
a = I

# ╔═╡ 44e5d727-392c-43bb-8dd1-f51dac984399
b = Zabck_

# ╔═╡ 1e10695a-f26e-44cb-990c-eec08f5a7d2e
c = zeros(Complex, 3, 3)

# ╔═╡ 7e35a316-27a3-4f6d-b754-58973397074c
d = a

# ╔═╡ 8a766ae0-d047-45e5-aa2d-11c42e061855
A = inv(a)

# ╔═╡ c651eea8-c6c0-40ed-b801-46d66ffb2a09
B = inv(a) * b

# ╔═╡ 08a8f081-d792-4a11-8d31-6bad523865dc
md"
## Definição de funções úteis

- p(m, a)
- dv(v)"

# ╔═╡ eb866aeb-1ced-48fd-b1ef-6a90ce813386
p(m, a) = m * cis(deg2rad(a))

# ╔═╡ ff8488b4-2970-44a3-bbbd-ad3608ff1ccb
function dv(v)
	for i in v
		m = abs(i)
		a = rad2deg(angle(i))
		@printf "%.2f ∠%.2fº\n" m a
	end
end

# ╔═╡ 9c28ac4e-094d-4f88-994d-59c469a253b9
md"""
## Definição das tensões do sistema e da carga conectada à Linha
"""

# ╔═╡ cc863d5c-32f1-4fc6-a074-feea2aa85c1b
v1l = 115e3

# ╔═╡ 7f4dec3a-bc68-40f2-97e1-4ff94958d791
v1f = v1l / √3

# ╔═╡ 1f6e74cf-01db-4245-9b52-4b4a90480636
v2l = 12.47e3

# ╔═╡ 9d013a7d-745e-4167-abd3-808ca24ba683
v2f = v2l / √3

# ╔═╡ 33d47c45-865d-44e5-b29e-9295af799b1a
Vabc1_ = [p(v1f, -30.0); p(v1f, -150.0); p(v1f, 90.0);;]

# ╔═╡ 06500d53-7961-4ab3-8671-4bf42cb0658c
# dv(Vabc1_)

# ╔═╡ ba8480a7-7e19-492d-bd77-6f19a080a91e
Sabc_ = [p(1.7e6, acosd(0.9)); p(1.2e6, acosd(0.85)); p(1.5e6, acosd(0.95));;]

# ╔═╡ f5889de3-685a-4c91-b109-380375fa7555
# dv(Sabc_)

# ╔═╡ 7327064c-e290-4f53-92df-d1001d4179a9


# ╔═╡ e3ecda18-2bf1-44ba-817b-d272cc1c75d4
md"""
## Matrizes para conversões de tensão

LN voltages:

$\mathbf{v}:=\left[\begin{array}{c}V_a \\ V_b \\ V_c\end{array}\right] \quad$

LL voltages:

$\quad \tilde{\mathbf{v}}:=\left[\begin{array}{c}V_{a b} \\ V_{b c} \\ V_{c a}\end{array}\right]$

LN-to-LL conversion:

$\quad \tilde{\mathbf{v}}=\mathbf{D}_f \mathbf{v}, \quad \mathbf{D}_f:=\left[\begin{array}{ccc}1 & -1 & 0 \\ 0 & 1 & -1 \\ -1 & 0 & 1\end{array}\right] \quad \begin{aligned} & \text { singular } \\ & \text { matrix! }\end{aligned}$

LL voltages are zero-sum:

$\quad \mathbf{1}^{\top} \tilde{\mathbf{v}}=V_{a b}+V_{b c}+V_{c a}=0 \quad \begin{gathered}\text { even for unbalanced } \\ \text { conditions }\end{gathered}$

Given LL voltages, recover equivalent LN voltages:

$\quad \mathbf{v}=\mathbf{W} \tilde{\mathbf{v}}, \quad \mathbf{W}:=\frac{1}{3}\left[\begin{array}{lll}2 & 1 & 0 \\ 0 & 2 & 1 \\ 1 & 0 & 2\end{array}\right]$

"""

# ╔═╡ 3ab2de54-50f9-44e0-aec3-045f3501cc08
md"""
## Matrizes para conversão de corretes

Line currents:

$\quad \mathbf{i}:=\left[\begin{array}{c}I_a \\ I_b \\ I_c\end{array}\right]$ 

Phase currents:

$\quad \tilde{\mathbf{i}}:=\left[\begin{array}{c}I_{a b} \\ I_{b c} \\ I_{c a}\end{array}\right]$

Phase to line conversion

$\mathbf{i}=\left[\begin{array}{c}
I_{a b}-I_{c a} \\
I_{b c}-I_{a b} \\
I_{c a}-I_{b c}
\end{array}\right]=\left[\begin{array}{ccc}
1 & 0 & -1 \\
-1 & 1 & 0 \\
0 & -1 & 1
\end{array}\right] \tilde{\mathbf{i}}=\mathbf{D}_f^{\top} \tilde{\mathbf{i}}$

If line currents exit triangle (delta source), vector $\mathbf{i}$ gets negative sign or $\tilde{\mathbf{i}}:=\left[\begin{array}{c}I_{b a} \\ I_{c b} \\ I_{a c}\end{array}\right]$;

Singularity (shift-invariance) can be waived by fixing the sum of delta currents;

Given line currents, recover equivalent delta currents

$\tilde{\mathbf{i}}=\mathbf{W}^{\top} \mathbf{i}$


"""

# ╔═╡ fd79f9b4-e466-470f-ae50-37a34e709ae0
Df = [1 -1 0; 
	  0 1 -1;
	  -1 0 1]

# ╔═╡ 9d697f89-7557-489c-9aca-65ed18a63720
W = 1/3 * [2 1 0; 
	       0 2 1;
	       1 0 2]

# ╔═╡ f6330451-b95f-4272-b727-675dbd7d17ec
Av = [0 1 0; 
	  0 0 1;
	  1 0 0]

# ╔═╡ f52c0958-14a6-469c-a836-0338b32ffe4d
md"""
## Definição das matrizes características do trafo 3ph $\Delta-Y$ step-down
"""

# ╔═╡ c8c5a634-837d-4a6a-ab85-232db346feff
nt_ = v1l / v2f # relação de transformação do trafo 1ph equivalente

# ╔═╡ 9e4f0e6c-990a-417f-a0d5-33814f269380
Ztpu = p(0.085, 85) # impedância em pu do trafo

# ╔═╡ 65f1c91f-e789-4b8c-b173-65049b7cb02d
Stnom = 5e6 # Potência de base do trafo

# ╔═╡ 97b60e7b-7404-473a-8200-298721ce1bf6
Ztbase = v2l^2 / Stnom # Impedância de base do trafo

# ╔═╡ 2cc24336-e677-4703-89cf-8179337fcf77
Zt = Ztpu * Ztbase # impedância do trafo

# ╔═╡ 1e707057-3754-4a62-af60-9cb2ba640b5a


# ╔═╡ 1d9114ab-7aa3-4449-bd46-e968946620a7
md"""
Matrizes de (a, b, c, d, A, B) do trafo, em função das matrizes Df, W e Av:

$\begin{aligned}
a_t & = - nt \cdot W \cdot A_v\\
b_t & = -nt \cdot Z_t \cdot W \cdot A_v\\
c_t & = \mathbf{0}\\
d_t & = \frac{1}{n_t} \cdot D_f^{'}\\
A_t & = d_t\\
B_t & = Z_t \cdot I
\end{aligned}$
"""

# ╔═╡ c78f2d08-378e-46d6-b2f7-20fb761c5548
ct = zeros(3, 3)

# ╔═╡ ff248426-dc21-4ea0-8546-35eba1c25680
Bt = Zt * I(3)

# ╔═╡ 4cab6751-d051-458b-a76a-04f713eb8320
# dv(Vabc1_)

# ╔═╡ df17a1b2-946f-4a7c-a2cf-d1377e7a94dd
md"## Execução do Fluxo de carga"

# ╔═╡ b6ac79e1-65ea-4d2f-a921-1f1c754f3125
md"""
Algoritmo de fluxo de carga de varredura direta inversa com inclusão do regulador de tensão:

```
        (1)                 (2)              (3)
 	     |                   |                |
  (~)---- -------( )( )------ ---------------- ---->
	     |                   |                |
   Fonte	Transformador        Linha de           Carga
		    D-Y                  Distribuição
```

"""

# ╔═╡ 6225cce3-a0ba-4d62-9c41-3e4dcc85f474
md"""
Execução de Fluxo de Carga utilizando método de varredura Direta-Inversa

> Varredura Inversa $\Longleftarrow$:
> 
> $I_{abc}^{(n)} = \mathbf{c}~V_{abc}^{(m)} + \mathbf{d}~I_{abc}^{(m)}$

> Varredura Direta $\Longrightarrow$:
> 
> $V_{abc}^{(m)} = \mathbf{A}~V_{abc}^{(n)} - \mathbf{B}~I_{abc}^{(m)}$
> 
> Em que: $\mathbf{A} = \mathbf{a}^{-1}$ e $\mathbf{B} = \mathbf{a}^{-1} \mathbf{b}$

.

"""

# ╔═╡ 73a5a3df-84da-48f7-bd79-861221de3790
md"""
## Execução de resultados
"""

# ╔═╡ c85e95e2-86ce-492a-be61-db6d23aba4a7


# ╔═╡ 5c07bada-cd49-4a9a-bf14-bc8670a7df2a
md"""
### Tensões
"""

# ╔═╡ 2230af12-b138-4df7-93ee-a049650e0d64
# dv(Vabc2k)

# ╔═╡ ec4e0061-979f-4e90-9536-bcc52c0dbc2e
# dv(Vabc3k)

# ╔═╡ c0aed8d3-cb33-4050-85bd-a015bb69934d
md"""
As tensões na carga para a base de 120V:
"""

# ╔═╡ 745864eb-1d7e-4a6b-8923-dc01e2e875d8
md"""
Alteração dos taps do transformador:
"""

# ╔═╡ 6f3feaf1-63ea-4ee9-83a7-67c525e2e7f5
@bind nt Slider( 0.7*nt_:0.05*nt_:1.3*nt_, default=nt_ )

# ╔═╡ d0c13908-dc05-4cc2-80fb-7d4dc7c212ca
at = -nt * (W * Av)

# ╔═╡ fcc0c8f8-3833-4967-884e-fc4398748eb7
bt = -nt * Zt * (W * Av)

# ╔═╡ 292b0b71-16c5-4f4e-930b-5b09d6c19d7a
dt = 1/nt * Df'

# ╔═╡ c6f4209f-0b89-465d-a64f-01d691e209ee
At = 1 / nt *  Df'

# ╔═╡ f7cd27e8-94d2-4738-a80b-dbd707adc027
nt

# ╔═╡ c95807d0-5b87-4bed-8ec0-8b90fade70be


# ╔═╡ 0b92fb96-7de1-4052-8229-e8673f0dbfc3
md"""
### Correntes
"""

# ╔═╡ 1a2c6945-e2d5-4a70-b80c-02039b13fd9a
# dv(Iabc2k)

# ╔═╡ a06b318b-546b-472b-b58e-0f1ccd9d2566
# dv(Iabc1k)

# ╔═╡ 386dd884-7a43-42e8-9aaa-4b1e5a89dbe2


# ╔═╡ 3113877c-def3-4386-9161-24660138d2a1
md"""
### Alteração na condição de balanceamento da carga
"""

# ╔═╡ 761755ac-9bed-4223-84d4-4ca69e360776
@bind aft Slider(0.0:0.5:5.0, default=1.0)

# ╔═╡ 1340dedb-2bf7-4964-aaa5-7c9c527b5f58
aft

# ╔═╡ 94f92cb2-2ba4-49a2-b464-f244236d8b41
Sabc = Sabc_ .* [aft; 1.0; 1.0/aft]

# ╔═╡ 960ab441-7b03-44b2-a49e-9a57d9dc7dc0
begin
	Vabc1k = Vabc1_
	Vabc2k = At * Vabc1_
	Vabc3k_ = Vabc3k = Vabc2k
	
	nk = 0
	
	while nk < 100
		global Vabc1k, Iabc1k, Vabc2k, Iabc2k, Vabc3k, Iabc3k, Vabc3k_, nk
		nk = nk + 1
		Iabc3k = conj.(Sabc ./ Vabc3k)

		##################
		# Backward Sweep #
		#################
		# Vabc2k = a * Vabc3k + b * Iabc3k
		Iabc2k = c * Vabc3k + d * Iabc3k

		# Vabc1k = at * Vabc2k + bt * Iabc2k
		Iabc1k = ct * Vabc2k + dt * Iabc2k

		#################
		# Forward Sweep #
		################
		Vabc1k = Vabc1_
		Vabc2k = At * Vabc1k - Bt * Iabc2k
		Vabc3k = A * Vabc2k - B * Iabc3k

		# verificação da condição de parada
		error = maximum(abs.((Vabc3k - Vabc3k_) ./ Vabc3k_))
		# @printf "Iteration %i absolute error = %.4f\n" nk error
		if error < 1e-6
			break
		else
			Vabc3k_ = Vabc3k
		end
		
	end
end

# ╔═╡ 3fb9ff76-f729-46ac-81dd-d95b254cb41d
# dv(Vabc3k / (7200.0/120.0))
abs.(Vabc3k / (7200.0/120.0))

# ╔═╡ f53ef149-1c11-41ee-87ce-88e6606b9afd
abs(sum(Iabc2k))

# ╔═╡ 5967fe9c-e4b8-40f3-84e6-cef59042ad14
abs(sum(Iabc1k))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "123ea9c1f46e4a6d19d315365279dde6c5cd209e"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─deb8bd3d-cffd-49af-8d16-35122a923ebc
# ╟─3dcaa6fe-080a-4fb1-b18c-871249e513d6
# ╠═3add8c37-e20b-48ca-b6b7-db6d3e38e6f9
# ╠═d36ab22c-2897-4dad-a693-c70eaa1f81dc
# ╟─2ff6e94d-872a-471a-9ce3-dfe2a491f31e
# ╠═16bcfca2-557b-4512-a4b5-7e05c98382cd
# ╠═623e3772-f53a-4765-a365-1117aa5110ae
# ╠═960141d1-7b14-4e05-8e6a-aa37484d2b8b
# ╟─664eb69b-9fc8-4c98-a899-40aac945d3de
# ╟─d4c98edb-74ff-46ab-b8b0-d7acb115cfc7
# ╠═d6b1dc9c-c3a9-46ea-ac64-937f31f35e4e
# ╠═248d7de5-6b58-495b-a788-2450706c5b9f
# ╟─478ddb02-083d-45a0-9aad-6459235a878a
# ╠═571affb9-0c8f-4541-a4b1-3bea3e45b036
# ╟─542462a7-8c53-422d-96ba-d4ce6db0deda
# ╠═2de5a128-85b8-4111-9dd5-5c6450d7ad1b
# ╟─e45f3e9f-71c6-4a63-90a0-ee635dec4a43
# ╠═a51756e4-d14f-430f-b18b-5beda9f321cf
# ╠═44e5d727-392c-43bb-8dd1-f51dac984399
# ╠═1e10695a-f26e-44cb-990c-eec08f5a7d2e
# ╠═7e35a316-27a3-4f6d-b754-58973397074c
# ╠═8a766ae0-d047-45e5-aa2d-11c42e061855
# ╠═c651eea8-c6c0-40ed-b801-46d66ffb2a09
# ╟─08a8f081-d792-4a11-8d31-6bad523865dc
# ╠═eb866aeb-1ced-48fd-b1ef-6a90ce813386
# ╠═ff8488b4-2970-44a3-bbbd-ad3608ff1ccb
# ╟─9c28ac4e-094d-4f88-994d-59c469a253b9
# ╠═cc863d5c-32f1-4fc6-a074-feea2aa85c1b
# ╠═7f4dec3a-bc68-40f2-97e1-4ff94958d791
# ╠═1f6e74cf-01db-4245-9b52-4b4a90480636
# ╠═9d013a7d-745e-4167-abd3-808ca24ba683
# ╠═33d47c45-865d-44e5-b29e-9295af799b1a
# ╠═06500d53-7961-4ab3-8671-4bf42cb0658c
# ╠═ba8480a7-7e19-492d-bd77-6f19a080a91e
# ╠═f5889de3-685a-4c91-b109-380375fa7555
# ╟─7327064c-e290-4f53-92df-d1001d4179a9
# ╟─e3ecda18-2bf1-44ba-817b-d272cc1c75d4
# ╟─3ab2de54-50f9-44e0-aec3-045f3501cc08
# ╠═fd79f9b4-e466-470f-ae50-37a34e709ae0
# ╠═9d697f89-7557-489c-9aca-65ed18a63720
# ╠═f6330451-b95f-4272-b727-675dbd7d17ec
# ╟─f52c0958-14a6-469c-a836-0338b32ffe4d
# ╠═c8c5a634-837d-4a6a-ab85-232db346feff
# ╠═9e4f0e6c-990a-417f-a0d5-33814f269380
# ╠═65f1c91f-e789-4b8c-b173-65049b7cb02d
# ╠═97b60e7b-7404-473a-8200-298721ce1bf6
# ╠═2cc24336-e677-4703-89cf-8179337fcf77
# ╟─1e707057-3754-4a62-af60-9cb2ba640b5a
# ╟─1d9114ab-7aa3-4449-bd46-e968946620a7
# ╠═d0c13908-dc05-4cc2-80fb-7d4dc7c212ca
# ╠═fcc0c8f8-3833-4967-884e-fc4398748eb7
# ╠═c78f2d08-378e-46d6-b2f7-20fb761c5548
# ╠═292b0b71-16c5-4f4e-930b-5b09d6c19d7a
# ╠═c6f4209f-0b89-465d-a64f-01d691e209ee
# ╠═ff248426-dc21-4ea0-8546-35eba1c25680
# ╠═4cab6751-d051-458b-a76a-04f713eb8320
# ╟─df17a1b2-946f-4a7c-a2cf-d1377e7a94dd
# ╟─b6ac79e1-65ea-4d2f-a921-1f1c754f3125
# ╟─6225cce3-a0ba-4d62-9c41-3e4dcc85f474
# ╠═960ab441-7b03-44b2-a49e-9a57d9dc7dc0
# ╟─73a5a3df-84da-48f7-bd79-861221de3790
# ╟─c85e95e2-86ce-492a-be61-db6d23aba4a7
# ╟─5c07bada-cd49-4a9a-bf14-bc8670a7df2a
# ╠═2230af12-b138-4df7-93ee-a049650e0d64
# ╠═ec4e0061-979f-4e90-9536-bcc52c0dbc2e
# ╟─c0aed8d3-cb33-4050-85bd-a015bb69934d
# ╠═3fb9ff76-f729-46ac-81dd-d95b254cb41d
# ╟─745864eb-1d7e-4a6b-8923-dc01e2e875d8
# ╠═6f3feaf1-63ea-4ee9-83a7-67c525e2e7f5
# ╠═f7cd27e8-94d2-4738-a80b-dbd707adc027
# ╟─c95807d0-5b87-4bed-8ec0-8b90fade70be
# ╟─0b92fb96-7de1-4052-8229-e8673f0dbfc3
# ╠═1a2c6945-e2d5-4a70-b80c-02039b13fd9a
# ╠═a06b318b-546b-472b-b58e-0f1ccd9d2566
# ╠═f53ef149-1c11-41ee-87ce-88e6606b9afd
# ╠═5967fe9c-e4b8-40f3-84e6-cef59042ad14
# ╟─386dd884-7a43-42e8-9aaa-4b1e5a89dbe2
# ╟─3113877c-def3-4386-9161-24660138d2a1
# ╠═761755ac-9bed-4223-84d4-4ca69e360776
# ╠═1340dedb-2bf7-4964-aaa5-7c9c527b5f58
# ╠═94f92cb2-2ba4-49a2-b464-f244236d8b41
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
