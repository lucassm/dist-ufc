### A Pluto.jl notebook ###
# v0.19.42

#> [frontmatter]
#> title = "Voltage Regulation"
#> date = "2024-08-12"
#> tags = ["power distribution", "power flow", "voltage regulator"]
#> description = "Power Distribution"
#> 
#>     [[frontmatter.author]]
#>     name = "Lucas Melo"
#>     url = "https://lucassm.pro"

using Markdown
using InteractiveUtils

# ╔═╡ 2ac9ef60-48af-4469-823f-19f39b02fd5a
begin
	using LinearAlgebra
	using Statistics
	using Printf
	using PlutoUI
end

# ╔═╡ f30766d2-60c8-4609-a4c1-08ae518cdea2
md"""
# Matrizes de Impedância e Fluxo de Potência com e sem Regulador de Tensão  
"""

# ╔═╡ 718c1b14-3a97-4577-a91d-5df044773339
PlutoUI.TableOfContents(include_definitions=false)

# ╔═╡ 14b1d201-ca7d-4d70-aae1-ea33d6b9b008
md"""
## Montagem da Matriz de Impedâncias Primitiva
"""

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

# ╔═╡ a08bf787-5af0-41ba-9416-dec1c61869a5
md"""
## Montagem da Matriz de Impedâncias de Fase (Redução de Kron)
"""

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

# ╔═╡ cfd29b88-f19a-49f5-b4f8-3a5d58a152cb
md"""
## Montagem da Matriz de Impedâncias de Sequência
"""

# ╔═╡ b7558e5a-7532-485d-aef5-6ab66f07fc19
begin
	zs = (Zabck[1, 1] + Zabck[2, 2] + Zabck[3, 3]) / 3.0
	zm 	= (Zabck[1, 2] + Zabck[1, 3] + Zabck[2, 1]) / 3.0
	Zabcks = [zs zm zm; zm zs zm; zm zm zs]
end

# ╔═╡ bb455147-84c9-4d61-95e7-e1c3dad00ec3
md"""
## Cálculo das Matrizes características da Linha
"""

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
md"## Execução do Fluxo de carga usando Matriz Reduzida de Kron"

# ╔═╡ 1cb57bcc-1714-4712-aa6f-adfd9fde98ad
md"Definição de duas funções úteis:
- p(m, a)
- dv(v)"

# ╔═╡ b3a75a45-7aa6-4b76-9fe0-95bf3236e7fc
p(m, a) = m * cis(deg2rad(a))

# ╔═╡ 40169f3f-55e2-4a95-96ff-be550d9fc815
as_ = p(1.0, 120.0)

# ╔═╡ 9005d12e-67f8-458e-9a68-9c9d9cc0786b
As = [1.0 1.0 1.0; 1.0 as_^2 as_; 1.0 as_ as_^2]

# ╔═╡ 2e4f11f8-c22f-4153-b71f-ab8b80a4a421
Z012 = inv(As) * Zabck *  As

# ╔═╡ 4618fbce-6fa2-4dcd-bd8b-b83598bc4706
Z012s = inv(As) * Zabcks *  As

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

# ╔═╡ fea16ff9-4b56-4f32-9fa0-04274d0ae760
Sabc = [p(2.5e6, acosd(0.9)); p(2.0e6, acosd(0.85)); p(1.5e6, acosd(0.95))]

# ╔═╡ 8f68a404-c44e-4b62-90f5-2cb63035c9e2
dv(Sabc)

# ╔═╡ 028e07c7-cce4-4bfa-b79b-5da3683f77e6
md"""
Execução de Fluxo de Carga utilizando método de varredura Direta-Inversa

Varredura Inversa (<---):

$I_{abc}^{(n)} = \mathbf{c}~V_{abc}^{(m)} + \mathbf{d}~I_{abc}^{(m)}$

Varredura Direta (--->):

$V_{abc}^{(m)} = \mathbf{A}~V_{abc}^{(n)} - \mathbf{B}~I_{abc}^{(m)}$

Em que:

$\mathbf{A} = \mathbf{a}^{-1}$

e

$\mathbf{B} = \mathbf{a}^{-1} \mathbf{b}$
"""

# ╔═╡ 3dcdb752-c54b-4a6a-85a8-c62e6324c064
md"""
Algoritmo de fluxo de carga de varredura direta inversa:

```
        (1)                 (2)              
 	     |                   |               
 (~)----- ------------------- ---->
	     |                   |
Fonte		Linha de           Carga
	    	Distribuição
```

"""

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

# ╔═╡ e5e99775-94e2-4975-999a-ff56d6e191d1
md"""
## Fluxo de carga com regulador de tensão na fonte
"""

# ╔═╡ fb57fc9a-c025-49df-8d62-bbae8f853630
md"""
Relação de transformação do TP:
"""

# ╔═╡ 72c46f0e-92b0-47ae-93da-281993bd44ae
NPT = 7200.0 / 120

# ╔═╡ 897dcaf1-8f3f-4d0e-93b0-9ecf8a31343c
md"""
Relação de transformação do TC:
"""

# ╔═╡ f86a074f-49bf-4630-a8e0-db3783017d26
begin
	CTp = 600.0
	CTs = 5.0
	CT = CTp/CTs
end

# ╔═╡ 4e5d223b-f19c-4654-a040-3812c7f8fdb9
md"""
As tensões encontradas no fluxo de carga sem regulador de tensão na base de 120V serão:
"""

# ╔═╡ ec6923da-5a42-40ea-8870-e8a07a2d3921
dv(Vabc2/NPT)

# ╔═╡ 7694831f-08ab-45cc-816c-0a812ff2ca50
VL = 120.0 # voltage level

# ╔═╡ a5f1a705-f183-4759-910c-77276d6e464d
BW = 2.0 # bandwidth

# ╔═╡ 9ab0abc6-d735-4005-b7f4-c8ce85ff0de7
V_target = VL - BW/2.0

# ╔═╡ 36dc7f5f-f307-49d0-b644-3df5364f9bb9
md"""
Impedância da linha em $\Omega$:
"""

# ╔═╡ 9369b076-dc8c-434d-9446-c69ce400fe44
Zline = (Vabc1 - Vabc2) ./ Iabc2

# ╔═╡ 8da13bd0-bd7b-41e1-9d90-7141d0f05ba7
dv(Zline)

# ╔═╡ 6bd4fff9-4f2c-47b9-b2c2-82c0df370b8a
md"""
Impedância média:
"""

# ╔═╡ b4a9d207-4300-40fa-a918-97a67ec02114
Zline_avg = mean(Zline)

# ╔═╡ aff083c7-807c-4035-902b-315a6736fb84
md"""
Cálculo do ajuste do **compensador de linha** em *Volts*:
"""

# ╔═╡ a309467e-e72f-454c-92bc-8f5d8d0f0746
ZLC_volts = Zline_avg * CTp / NPT

# ╔═╡ b747f862-fa6e-446e-bf62-050d9e8055c4
taps = round.((V_target .- abs.(Vabc2)/NPT) / 0.75)

# ╔═╡ 3cd1315c-b813-4d13-b7ac-c3eddf5eda5e
md"""
!!! note "Por que a divisão por 0,75?"

	Se considerarmos a base de 120V e um ajuste de +/- 10% na tensão em 16 passos, teremos  a seguinte razão:

	$\frac{0,1 \times 120}{16} = \frac{12}{16} = 0,75$
"""

# ╔═╡ cd721743-8396-4e3d-8d10-2e116bf576fd
md"""
Definição das matrizes características do regulador de tensão:
"""

# ╔═╡ 6affdbff-18d4-4de3-9533-d8aa505599ec
areg = I - 0.00625 * diagm(vec(taps))

# ╔═╡ 45d8a352-8a86-452f-b62a-1813a5276ccb
breg = creg = zeros(3, 3)

# ╔═╡ b164859e-e159-42c3-802b-c0dbded5e01d
dreg = Areg = inv(areg)

# ╔═╡ 81a52d73-a089-45f8-8f0b-2937669d199d
Breg = zeros(3, 3)

# ╔═╡ 3dd6a1a6-44c9-43d4-8695-31fe99dca660
dv(Vabc1_)

# ╔═╡ 7d570ef7-7f94-4b10-a4a6-14e884246af7
md"""
Algoritmo de fluxo de carga de varredura direta inversa com inclusão do regulador de tensão:

```
        (1)                 (2)              (3)
 	     |                   |                |
  (~)---- --------()(/)------ ---------------- ---->
	     |                   |                |
   Fonte		Regulador        Linha de           Carga
		    	de Tensão        Distribuição
```

"""

# ╔═╡ 476c8ad7-b74a-4fcc-968d-261a309edc9f
begin
	Vabc3k = Vabc1_
	Vabc2ar = Vabc1_
	nk = 0
	
	while nk < 100
		global Vabc1k, Iabc1k, Vabc2ar, Iabc2ar, Vabc3k, Iabc3k, nk
		nk = nk + 1
		Iabc3k = conj.(Sabc ./ Vabc3k)

		# Backward Sweep
		Vabc2ar = a * Vabc3k + b * Iabc3k
		Iabc2ar = c * Vabc3k + d * Iabc3k

		Vabc1k = areg * Vabc2ar + breg * Iabc2ar
		Iabc1k = creg * Vabc2ar + dreg * Iabc2ar
		

		error = maximum(abs.((Vabc1k - Vabc1_) ./ Vabc1_))
		@printf "Iteration %i absolute error = %.4f\n" nk error
		if error < 1e-5
			break
		end

		# Forward Sweep
		Vabc1k = Vabc1_
		Vabc2ar = Areg * Vabc1k - Breg * Iabc2ar
		Vabc3k = A * Vabc2ar - B * Iabc3k
	end
end

# ╔═╡ c4d6922d-9674-479f-99a7-eab31d606a3e
dv(Vabc3k / NPT)

# ╔═╡ 2c45b016-6a49-4095-9362-854626fc27ca
dv(Vabc1_)

# ╔═╡ 1f46acc9-b0e3-43e2-860c-6e0505665c6f
dv(Iabc3k)

# ╔═╡ 5ebfff46-55ca-4114-9eeb-12134a2e9a9f
dv(Iabc1k)

# ╔═╡ e5f73b84-521b-46dc-89e1-b8d65fd9a443
md"""
!!! danger "Importante!" 
	Note que em nenhum momento o valor da **impedância do compensador de linha** foi utilizado neste exemplo!

	Esse fato gera uma **imprecisão nos cálculos**, conforme demonstrado abaixo.
"""

# ╔═╡ 74ce4fc7-379d-4f8a-b160-453c38fbfce3
Zc = ZLC_volts / CTs

# ╔═╡ 60ca5065-3ac1-4127-bb94-de7f7bd03519
# Para os taps do regulador inicialmente na posição 0
# areg = I
V_reg_out = I * Vabc1_

# ╔═╡ f09d4b2a-5fce-4fae-9d77-3b0790ce3421
dv(V_reg_out)

# ╔═╡ ec97f41d-442f-4d35-81c8-fbf997bcb778
# Para os taps do regulador inicialmente na posição 0
# dreg = inv(I)
I_reg_out = inv(I) * Iabc2

# ╔═╡ 56781a42-b807-4910-bd5b-f73ccee7e504
dv(I_reg_out)

# ╔═╡ fe61f974-14bb-4927-95c4-420e126013f0
V_reg_s = V_reg_out / NPT

# ╔═╡ 3a3d836d-1934-4f16-9639-84f050c5264d
dv(V_reg_s)

# ╔═╡ 203bf3f3-dca5-4a1c-8cc6-207945ade0f6
I_reg_s = I_reg_out / CT

# ╔═╡ f042cc86-4e7c-47ae-b6d3-9bbc6015e49e
dv(I_reg_s)

# ╔═╡ 1ad1786b-b4c6-44c8-be16-f5847012f40b
Zcomp_matrix = Zc * I(3)

# ╔═╡ 9f3190fe-6371-4a84-80f8-a332daf3701b
V_relay = V_reg_s - Zcomp_matrix * I_reg_s

# ╔═╡ 2c821fd7-e3fe-4406-a83b-b1436c331d43
dv(V_relay)

# ╔═╡ 7094dc6c-fb8e-4f9a-93db-3bed08da858f
taps_ = round.((V_target .- abs.(V_relay)) / 0.75)

# ╔═╡ 3f3225a4-1377-4019-83e9-444c21d261dd
md"""
!!! warning "Importante!"
	Veja como esses valores de tap são diferentes dos valores utilizados acima para o cálculo do fluxo de carga!
"""

# ╔═╡ 1ec0be9a-cb64-4946-b0b8-a8ece0bbbfaa
# TODO: Executar fluxo de carga com os valores de tap
# do regulador de tensão presnetes na variável taps_

# ╔═╡ e51fd0a3-3e02-4fda-8862-93c9fea315e8
md"""
## Cálculo de fluxo de carga com análise dos taps um por vez
"""

# ╔═╡ f982a709-95ad-46f0-8e6f-635bab0b7b32
begin

	taps_h = [0; 0; 0]

	h = 0
	while h < 16
		global Vabc1_h, Vabc2ar_h, Vabc3_h, h, nh, taps_h
		h = h + 1
		
		areg_h = I - 0.00625 * diagm(vec(taps_h))
		breg_h = creg_h = zeros(3, 3)
		dreg_h = Areg_h = inv(areg_h)
		Breg_h = zeros(3, 3)
		
		Vabc3_h = Vabc1_
		Vabc2ar_h = Vabc1_
		Vabc1_h = Vabc1_
		
		nh = 0
		
		while nh < 15
			global Vabc1_h, Vabc2ar_h, Vabc3_h, Iabc1_h, Iabc2ar_h, Iabc3_h, nh
			nh = nh + 1
			Iabc3_h = conj.(Sabc ./ Vabc3_h)
	
			# Backward Sweep
			Vabc2ar_h = a * Vabc3_h + b * Iabc3_h
			Iabc2ar_h = c * Vabc3_h + d * Iabc3_h
	
			Vabc1_h = areg_h * Vabc2ar_h + breg_h * Iabc2ar_h
			Iabc1_h = creg_h * Vabc2ar_h + dreg_h * Iabc2ar_h
			
	
			error = maximum(abs.((Vabc1_h - Vabc1_) ./ Vabc1_))
			@printf "Iteration %i absolute error = %.4f\n" nh error
			if error < 1e-5
				break
			end
	
			# Forward Sweep
			Vabc1_h = Vabc1_
			Vabc2ar_h = Areg_h * Vabc1_h - Breg_h * Iabc2ar_h
			Vabc3_h = A * Vabc2ar_h - B * Iabc3_h
		end

		# Cálculo da tensão do relé de atuação do LDC no 
		# circuito de comando do Regulador de Tensão
		I_reg_out_h = Iabc2ar_h
		
		I_reg_s_h = I_reg_out_h / CT
		V_reg_s_h  = Vabc2ar_h / NPT
		
		V_relay_h = V_reg_s_h - Zcomp_matrix * I_reg_s_h
		
		if sum(abs.(V_relay_h) .< V_target) == 0.0
			println("Taps: ", taps_h)
			dv(V_relay_h)
			break
		else
			println("Taps: ", taps_h)
			taps_h = taps_h + (abs.(V_relay_h) .< V_target)
			dv(V_relay_h)
		end
	end
end

# ╔═╡ 526359b6-520c-494b-81da-ca6e1ba5e1e9
taps_h

# ╔═╡ 0031b33e-8c15-4bc9-bd82-5495b50c068f
md"""
!!! warning "Tensão na carga após a atuação do regulador"
	Para o ajuste implementado, em decorrância da situação de carga elevada na fase A, a tensão nesta fase ficará fora da faixa de ajuste, conforme exibido abaixo:
"""

# ╔═╡ 7d7178f4-394e-441b-91c6-909163bb9d75
dv(Vabc3_h/NPT)

# ╔═╡ 436ea71b-25ba-4b4f-8d2b-e3bf68a7b257
md"""
!!! tip "Dica do autor para ajuste do regulador de tensão"
	 Because the three line currents are all different, it means the heavily loaded phase (a) voltage will not represent what is actually happening on the system. Once again, this is a problem that occurs because of the unbalanced loading.

	One way to raise the load voltages is: **to specify a higher voltage level by increasing the voltage level to 122 V**.
"""

# ╔═╡ 230f84ce-4dbc-4655-9f62-e56e763c33c1
begin
	struct TwoColumn{L, R}
	    left::L
	    right::R
	end
	
	function Base.show(io, mime::MIME"text/html", tc::TwoColumn)
	    write(io, """<div style="display: flex;"><div style="flex: 50%;">""")
	    show(io, mime, tc.left)
	    write(io, """</div><div style="flex: 50%;">""")
	    show(io, mime, tc.right)
	    write(io, """</div></div>""")
	end
end

# ╔═╡ d3a26f4e-0e6b-4aee-875e-0ec7bae04ef2


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

julia_version = "1.10.5"
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
version = "5.11.0+0"

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
# ╟─f30766d2-60c8-4609-a4c1-08ae518cdea2
# ╠═2ac9ef60-48af-4469-823f-19f39b02fd5a
# ╠═718c1b14-3a97-4577-a91d-5df044773339
# ╠═14b1d201-ca7d-4d70-aae1-ea33d6b9b008
# ╠═231f501e-4893-11ef-2712-61982ea6360c
# ╠═6ac960cb-768f-4a8d-8f00-cfee2e929710
# ╠═a12d9137-5191-4e8a-ab01-c77e7429ea22
# ╟─9664cdaf-afbc-491f-8740-dd13ecdd0c21
# ╠═dd221b04-6a9f-4329-84a7-992c90a7844e
# ╟─4aa82beb-a1c9-42f4-b18b-3a1947032eee
# ╠═1acc18f4-e082-4400-ad9e-111c0e9a47d4
# ╟─a08bf787-5af0-41ba-9416-dec1c61869a5
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
# ╟─cfd29b88-f19a-49f5-b4f8-3a5d58a152cb
# ╠═40169f3f-55e2-4a95-96ff-be550d9fc815
# ╠═9005d12e-67f8-458e-9a68-9c9d9cc0786b
# ╠═2e4f11f8-c22f-4153-b71f-ab8b80a4a421
# ╠═b7558e5a-7532-485d-aef5-6ab66f07fc19
# ╠═4618fbce-6fa2-4dcd-bd8b-b83598bc4706
# ╟─bb455147-84c9-4d61-95e7-e1c3dad00ec3
# ╠═a9839a94-86d2-4cf8-bca0-9927bcc2daf1
# ╠═7fd02a9f-573a-4e59-8cb1-2c20de856a22
# ╠═02168ec5-7999-4b19-99e8-3a349ffa7f4c
# ╠═ea84a6b5-6aa9-4aa5-9163-9914773d63e2
# ╠═6e182e50-be00-448b-92dd-07eee30ad05e
# ╠═e5930e39-2948-45cf-a0cc-453d51cb7c6c
# ╠═9abb5ad8-540f-48a8-a16d-7a7befaea988
# ╟─1cb57bcc-1714-4712-aa6f-adfd9fde98ad
# ╠═b3a75a45-7aa6-4b76-9fe0-95bf3236e7fc
# ╠═99f99038-171a-49c4-bdf0-149f21652e95
# ╠═2fd2761c-7252-456b-9133-42d15b62fc57
# ╠═f681de1a-0e63-4736-8a47-24c74018b5d6
# ╠═da1ca8fe-cdcd-416d-9268-02e4c548261f
# ╠═fea16ff9-4b56-4f32-9fa0-04274d0ae760
# ╠═8f68a404-c44e-4b62-90f5-2cb63035c9e2
# ╟─028e07c7-cce4-4bfa-b79b-5da3683f77e6
# ╟─3dcdb752-c54b-4a6a-85a8-c62e6324c064
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
# ╟─01e6a603-63ce-4369-93f6-4d5008e4b499
# ╠═8b1b27d2-00ff-4bcd-8a72-ef67e01a3c76
# ╠═055b2396-b55e-46f4-a331-1aeb7cd44550
# ╟─51cdf7ed-c72d-4056-8cf2-65212b86bfab
# ╠═d2af89d7-3dd1-4c19-9659-4edb2e6cba03
# ╟─e5e99775-94e2-4975-999a-ff56d6e191d1
# ╟─fb57fc9a-c025-49df-8d62-bbae8f853630
# ╠═72c46f0e-92b0-47ae-93da-281993bd44ae
# ╟─897dcaf1-8f3f-4d0e-93b0-9ecf8a31343c
# ╠═f86a074f-49bf-4630-a8e0-db3783017d26
# ╟─4e5d223b-f19c-4654-a040-3812c7f8fdb9
# ╠═ec6923da-5a42-40ea-8870-e8a07a2d3921
# ╠═7694831f-08ab-45cc-816c-0a812ff2ca50
# ╠═a5f1a705-f183-4759-910c-77276d6e464d
# ╠═9ab0abc6-d735-4005-b7f4-c8ce85ff0de7
# ╟─36dc7f5f-f307-49d0-b644-3df5364f9bb9
# ╠═9369b076-dc8c-434d-9446-c69ce400fe44
# ╠═8da13bd0-bd7b-41e1-9d90-7141d0f05ba7
# ╟─6bd4fff9-4f2c-47b9-b2c2-82c0df370b8a
# ╠═b4a9d207-4300-40fa-a918-97a67ec02114
# ╟─aff083c7-807c-4035-902b-315a6736fb84
# ╠═a309467e-e72f-454c-92bc-8f5d8d0f0746
# ╠═b747f862-fa6e-446e-bf62-050d9e8055c4
# ╟─3cd1315c-b813-4d13-b7ac-c3eddf5eda5e
# ╟─cd721743-8396-4e3d-8d10-2e116bf576fd
# ╠═6affdbff-18d4-4de3-9533-d8aa505599ec
# ╠═45d8a352-8a86-452f-b62a-1813a5276ccb
# ╠═b164859e-e159-42c3-802b-c0dbded5e01d
# ╠═81a52d73-a089-45f8-8f0b-2937669d199d
# ╠═3dd6a1a6-44c9-43d4-8695-31fe99dca660
# ╟─7d570ef7-7f94-4b10-a4a6-14e884246af7
# ╠═476c8ad7-b74a-4fcc-968d-261a309edc9f
# ╠═c4d6922d-9674-479f-99a7-eab31d606a3e
# ╠═2c45b016-6a49-4095-9362-854626fc27ca
# ╠═1f46acc9-b0e3-43e2-860c-6e0505665c6f
# ╠═5ebfff46-55ca-4114-9eeb-12134a2e9a9f
# ╟─e5f73b84-521b-46dc-89e1-b8d65fd9a443
# ╠═74ce4fc7-379d-4f8a-b160-453c38fbfce3
# ╠═60ca5065-3ac1-4127-bb94-de7f7bd03519
# ╠═f09d4b2a-5fce-4fae-9d77-3b0790ce3421
# ╠═ec97f41d-442f-4d35-81c8-fbf997bcb778
# ╠═56781a42-b807-4910-bd5b-f73ccee7e504
# ╠═fe61f974-14bb-4927-95c4-420e126013f0
# ╠═3a3d836d-1934-4f16-9639-84f050c5264d
# ╠═203bf3f3-dca5-4a1c-8cc6-207945ade0f6
# ╠═f042cc86-4e7c-47ae-b6d3-9bbc6015e49e
# ╠═1ad1786b-b4c6-44c8-be16-f5847012f40b
# ╠═9f3190fe-6371-4a84-80f8-a332daf3701b
# ╠═2c821fd7-e3fe-4406-a83b-b1436c331d43
# ╠═7094dc6c-fb8e-4f9a-93db-3bed08da858f
# ╟─3f3225a4-1377-4019-83e9-444c21d261dd
# ╠═1ec0be9a-cb64-4946-b0b8-a8ece0bbbfaa
# ╟─e51fd0a3-3e02-4fda-8862-93c9fea315e8
# ╠═f982a709-95ad-46f0-8e6f-635bab0b7b32
# ╠═526359b6-520c-494b-81da-ca6e1ba5e1e9
# ╟─0031b33e-8c15-4bc9-bd82-5495b50c068f
# ╠═7d7178f4-394e-441b-91c6-909163bb9d75
# ╟─436ea71b-25ba-4b4f-8d2b-e3bf68a7b257
# ╟─230f84ce-4dbc-4655-9f62-e56e763c33c1
# ╠═d3a26f4e-0e6b-4aee-875e-0ec7bae04ef2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
