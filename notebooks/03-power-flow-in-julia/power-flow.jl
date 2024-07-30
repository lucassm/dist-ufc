### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 2ac9ef60-48af-4469-823f-19f39b02fd5a
begin
	using LinearAlgebra
	using Printf
end

# ╔═╡ 231f501e-4893-11ef-2712-61982ea6360c
D = [0.0244 2.5 7.0 5.6569;
	 2.5 0.0244 4.5 4.272;
	 7.0 4.5 0.0244 5.0;
	 5.6569 4.272 5.0 0.0081]

# ╔═╡ 6ac960cb-768f-4a8d-8f00-cfee2e929710
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

# ╔═╡ a12d9137-5191-4e8a-ab01-c77e7429ea22
Zabc

# ╔═╡ 9664cdaf-afbc-491f-8740-dd13ecdd0c21
md"""
Conversão da matriz de impedâncias primitiva de Ohms por milha para Ohms por metro:
"""

# ╔═╡ dd221b04-6a9f-4329-84a7-992c90a7844e
Zabc * 0.000621371

# ╔═╡ 4aa82beb-a1c9-42f4-b18b-3a1947032eee
md"""
Conversão da matriz de impedâncias primitiva de Ohms por milha para Ohms por quilômetro:
"""

# ╔═╡ 1acc18f4-e082-4400-ad9e-111c0e9a47d4
Zabc * 0.621371

# ╔═╡ d1e1d24b-bc45-4bfd-b895-02cd04ea95bd
md"""
Redução de Kron da matriz de impedâncias primitiva
"""

# ╔═╡ e874b95f-5ea7-469d-ac4e-879efff9a006
Zabck = Zabc[1:3, 1:3] - Zabc[1:3, 4:4] * Zabc[4:4, 1:3] / Zabc[4, 4] 

# ╔═╡ ba31541d-dd91-40dd-8750-cae3a475656c
md"""
Conversão da matriz de impedâncias de fase de Ohms por milha para Ohms por metro:
"""

# ╔═╡ c532c853-98d8-40b1-a5b8-2daf89f3f2d9
Zabck * 0.000621371

# ╔═╡ 4ca75512-80d1-4a85-ab4b-f1616063d922
md"""
Conversão da matriz de impedâncias de fase de Ohms por milha para Ohms por quilômetro:

"""

# ╔═╡ eb3d79f8-883e-4c5c-bf79-fca4517303c4
Zabck * 0.621371

# ╔═╡ 87e7126f-7070-42fe-8767-8f76c1974de3
l = 10e3 / 5280.0

# ╔═╡ efa57c44-cfc5-4b2d-a9a3-7bbf8a9f4636
md"""
Matriz de Impedâncias Primitiva em Ohms:
"""

# ╔═╡ 311225df-49e4-4456-a752-0039f33dbba9
Zabc_ = Zabc * l

# ╔═╡ 79f0dadc-da16-48c3-a79c-18f9c8da8fb8
md"""
Matriz de Impedâncias de Fase Reduzida de Kron em Ohms:
"""

# ╔═╡ 9d58d71a-28fc-4695-b445-2215083fccf4
Zabck_ = Zabck * l

# ╔═╡ a9839a94-86d2-4cf8-bca0-9927bcc2daf1
a = I

# ╔═╡ 7fd02a9f-573a-4e59-8cb1-2c20de856a22
b = Zabck_

# ╔═╡ 02168ec5-7999-4b19-99e8-3a349ffa7f4c
c = zeros(Complex, 3, 3)

# ╔═╡ ea84a6b5-6aa9-4aa5-9163-9914773d63e2
d = a

# ╔═╡ 6e182e50-be00-448b-92dd-07eee30ad05e
A = inv(a)

# ╔═╡ e5930e39-2948-45cf-a0cc-453d51cb7c6c
B = A * b

# ╔═╡ 9abb5ad8-540f-48a8-a16d-7a7befaea988
md"# Execução do Fluxo de carga usando Matriz Reduzida de Kron"

# ╔═╡ 1cb57bcc-1714-4712-aa6f-adfd9fde98ad
md"Definição de duas funções úteis:
- p(m, a)
- dv(v)"

# ╔═╡ b3a75a45-7aa6-4b76-9fe0-95bf3236e7fc
p(m, a) = m * cis(deg2rad(a))

# ╔═╡ 99f99038-171a-49c4-bdf0-149f21652e95
function dv(v)
	for i in v
		m = abs(i)
		a = rad2deg(angle(i))
		@printf "%.2f ∠%.2fº\n" m a
	end
end

# ╔═╡ 2fd2761c-7252-456b-9133-42d15b62fc57
va = 12.47e3 / √3

# ╔═╡ f681de1a-0e63-4736-8a47-24c74018b5d6
Vabc1_ = [p(va, 0.0); p(va, -120.0); p(va, 120.0);;]

# ╔═╡ da1ca8fe-cdcd-416d-9268-02e4c548261f
dv(Vabc1_)

# ╔═╡ 70ff5ef4-23ca-4452-bd24-a1f21219f1d4
s = 2.0e6

# ╔═╡ fea16ff9-4b56-4f32-9fa0-04274d0ae760
Sabc = fill(p(s, acosd(0.9)), (3, 1))

# ╔═╡ 8f68a404-c44e-4b62-90f5-2cb63035c9e2
dv(Sabc)

# ╔═╡ f10df600-28c5-45fd-8e0d-f6c2909cbb95
begin
	Vabc2 = Vabc1_
	n = 0
	while n < 10
		global Vabc1, Iabc2, Vabc2, n
		n = n + 1
		Iabc2 = conj.(Sabc ./ Vabc2)

		# Backward Sweep
		Vabc1 = a * Vabc2 + b * Iabc2
		Iabc1 = c * Vabc2 + d * Iabc2

		error = maximum(abs.((Vabc1 - Vabc1_) ./ Vabc1_))
		@printf "Iteration %i absolute error = %.4f\n" n error
		if error < 1e-5
			break
		end

		# Forward Sweep
		Vabc1 = Vabc1_
		Vabc2 = A * Vabc1 - B * Iabc2
	end
end

# ╔═╡ 66ee744a-7147-44ad-8a23-94f301ce245b
md"Tensões na carga:"

# ╔═╡ 17eb3c0f-e0a7-4a21-8407-9044faa8d770
dv(Vabc2)

# ╔═╡ 4c49c86d-dc8a-4003-aa5f-200404cf043a
md"Correntes na linha:"

# ╔═╡ 28bd48fa-0f93-4030-ad1d-48984c60ba6c
dv(Iabc2)

# ╔═╡ b9d50eb8-a941-4016-8a8a-b57f6a59e199
md"Tensões na fonte:"

# ╔═╡ 135976c8-e727-4f8b-9bb2-a8d8f060e9e2
dv(Vabc1)

# ╔═╡ ea2ef140-8fa8-4950-9f6a-2178a8e62034
md"Razão entre as tensões na carga e na fonte:

$V_{abc}^{(2)} / V_{abc}^{(1)}$"

# ╔═╡ c9bd5963-06f0-47ee-863f-2ed07e66f186
abs.(Vabc2) ./ abs.(Vabc1)

# ╔═╡ 4954f6c8-1550-496a-b6e9-e9938ebf99c4
md"Corrente de desequilíbrio na carga:"

# ╔═╡ 7cfb48f7-8b50-4fed-8a50-4e89f1843342
abs(sum(Iabc2))

# ╔═╡ 49c5b75c-33a4-4c0f-a970-eaf9200e3982
md"Corrente no condutor neutro:"

# ╔═╡ 924fefd3-93f4-45ea-b414-d6540bf1b077
begin
	tn = - Zabc[4:4, 1:3] / Zabc[4, 4]
	In2 = tn * Iabc2
end

# ╔═╡ 9c937f87-3fa3-488b-a516-d1207f28e326
dv(In2)

# ╔═╡ 01e6a603-63ce-4369-93f6-4d5008e4b499
md"Corrente de retorno pela terra:"

# ╔═╡ 8b1b27d2-00ff-4bcd-8a72-ef67e01a3c76
Ig2 = -[sum(Iabc2)] - In2

# ╔═╡ 055b2396-b55e-46f4-a331-1aeb7cd44550
dv(Ig2)

# ╔═╡ 51cdf7ed-c72d-4056-8cf2-65212b86bfab
md"Corrente de Desequilíbrio total:"

# ╔═╡ d2af89d7-3dd1-4c19-9659-4edb2e6cba03
dv(In2 + Ig2)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "8d48bd265fda7892a322f49802c215af0234a337"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╠═2ac9ef60-48af-4469-823f-19f39b02fd5a
# ╠═231f501e-4893-11ef-2712-61982ea6360c
# ╠═6ac960cb-768f-4a8d-8f00-cfee2e929710
# ╠═a12d9137-5191-4e8a-ab01-c77e7429ea22
# ╟─9664cdaf-afbc-491f-8740-dd13ecdd0c21
# ╠═dd221b04-6a9f-4329-84a7-992c90a7844e
# ╟─4aa82beb-a1c9-42f4-b18b-3a1947032eee
# ╠═1acc18f4-e082-4400-ad9e-111c0e9a47d4
# ╟─d1e1d24b-bc45-4bfd-b895-02cd04ea95bd
# ╠═e874b95f-5ea7-469d-ac4e-879efff9a006
# ╟─ba31541d-dd91-40dd-8750-cae3a475656c
# ╠═c532c853-98d8-40b1-a5b8-2daf89f3f2d9
# ╟─4ca75512-80d1-4a85-ab4b-f1616063d922
# ╠═eb3d79f8-883e-4c5c-bf79-fca4517303c4
# ╠═87e7126f-7070-42fe-8767-8f76c1974de3
# ╟─efa57c44-cfc5-4b2d-a9a3-7bbf8a9f4636
# ╠═311225df-49e4-4456-a752-0039f33dbba9
# ╟─79f0dadc-da16-48c3-a79c-18f9c8da8fb8
# ╠═9d58d71a-28fc-4695-b445-2215083fccf4
# ╠═a9839a94-86d2-4cf8-bca0-9927bcc2daf1
# ╠═7fd02a9f-573a-4e59-8cb1-2c20de856a22
# ╠═02168ec5-7999-4b19-99e8-3a349ffa7f4c
# ╠═ea84a6b5-6aa9-4aa5-9163-9914773d63e2
# ╠═6e182e50-be00-448b-92dd-07eee30ad05e
# ╠═e5930e39-2948-45cf-a0cc-453d51cb7c6c
# ╟─9abb5ad8-540f-48a8-a16d-7a7befaea988
# ╟─1cb57bcc-1714-4712-aa6f-adfd9fde98ad
# ╠═b3a75a45-7aa6-4b76-9fe0-95bf3236e7fc
# ╠═99f99038-171a-49c4-bdf0-149f21652e95
# ╠═2fd2761c-7252-456b-9133-42d15b62fc57
# ╠═f681de1a-0e63-4736-8a47-24c74018b5d6
# ╠═da1ca8fe-cdcd-416d-9268-02e4c548261f
# ╠═70ff5ef4-23ca-4452-bd24-a1f21219f1d4
# ╠═fea16ff9-4b56-4f32-9fa0-04274d0ae760
# ╠═8f68a404-c44e-4b62-90f5-2cb63035c9e2
# ╠═f10df600-28c5-45fd-8e0d-f6c2909cbb95
# ╟─66ee744a-7147-44ad-8a23-94f301ce245b
# ╠═17eb3c0f-e0a7-4a21-8407-9044faa8d770
# ╟─4c49c86d-dc8a-4003-aa5f-200404cf043a
# ╠═28bd48fa-0f93-4030-ad1d-48984c60ba6c
# ╟─b9d50eb8-a941-4016-8a8a-b57f6a59e199
# ╠═135976c8-e727-4f8b-9bb2-a8d8f060e9e2
# ╟─ea2ef140-8fa8-4950-9f6a-2178a8e62034
# ╠═c9bd5963-06f0-47ee-863f-2ed07e66f186
# ╟─4954f6c8-1550-496a-b6e9-e9938ebf99c4
# ╠═7cfb48f7-8b50-4fed-8a50-4e89f1843342
# ╟─49c5b75c-33a4-4c0f-a970-eaf9200e3982
# ╠═924fefd3-93f4-45ea-b414-d6540bf1b077
# ╠═9c937f87-3fa3-488b-a516-d1207f28e326
# ╠═01e6a603-63ce-4369-93f6-4d5008e4b499
# ╠═8b1b27d2-00ff-4bcd-8a72-ef67e01a3c76
# ╠═055b2396-b55e-46f4-a331-1aeb7cd44550
# ╠═51cdf7ed-c72d-4056-8cf2-65212b86bfab
# ╠═d2af89d7-3dd1-4c19-9659-4edb2e6cba03
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
