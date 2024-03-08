### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ ecc1218c-a408-4de6-8cb5-40d29613b9e7
# ╠═╡ show_logs = false
import Pkg; Pkg.activate("/home/iwsatlas1/henkes/l200/auto/")

# ╔═╡ 2f337776-68cd-4b38-9047-88742bfa1c8a
begin
	using PlutoUI
	import PlutoUI: combine
end

# ╔═╡ 493bdebf-7802-4938-8693-5e0648bd8a2b
# ╠═╡ show_logs = false
begin
	ENV["JULIA_DEBUG"] = Main # enable debug
	ENV["JULIA_CPU_TARGET"] = "generic" # enable AVX2
	# import Pkg; Pkg.instantiate(); Pkg.precompile() # load packages
	
	using LegendDataManagement, PropertyFunctions, TypedTables, PropDicts
	using Unitful, Formatting, LaTeXStrings
	using LegendHDF5IO, LegendDSP, LegendSpecFits
	using Distributed, ProgressMeter

	using LegendDataTypes
	using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked

	using Measures
	using StatsBase
	
	l200 = LegendData(:l200)
	@info "Loading Legend MetaData"
end

# ╔═╡ cf295ad2-13ca-4508-bf51-5cfc52229fbd
# ╠═╡ show_logs = false
begin
	using Plots
	plotlyjs()
end;

# ╔═╡ 98a808df-40ff-4eda-89ee-772147f7e42f
html"""<style>
main {
	max-width: 80%;
	margin-right: 0;
}
"""

# ╔═╡ f2cf0e3f-d544-4a0b-bcdf-d02bf9beb9d6
html"""
<style>
input[type*="range"] {
	width: 30%;
}
</style>
"""

# ╔═╡ d35a8320-ae1e-485d-8751-1a2b36f2b809
legend_logo = "https://legend-exp.org/typo3conf/ext/sitepackage/Resources/Public/Images/Logo/logo_legend_tag_next.svg";

# ╔═╡ 4933bfe7-2d34-4283-bd4f-2cb0006c5158
md"""$(Resource(legend_logo))
## Julia Analysis Software - Energy Calibration Inspector
This tool lets you investigate a certain set of waveforms for a specififc detector and how certain parameters impact the DSP.
"""

# ╔═╡ 24b387ff-5df5-45f9-8857-830bc6580857
md""" ### Parameter settings
"""

# ╔═╡ 39b468d4-673c-4834-b042-d80cd9e0dfe7
md"Select Period"

# ╔═╡ 2c4dd859-0aa6-4cd2-94b3-7a171368065b
begin
	periods = search_disk(DataPeriod, l200.tier[:raw, :cal])
	@bind period Select(periods, default=DataPeriod(3))
end

# ╔═╡ 63b6be12-3b55-450a-8a2a-3d410f0dcb6c
md"Select Run"

# ╔═╡ 36232152-9811-4520-ba65-0f855e4ebbe4
begin
	runs = search_disk(DataRun, l200.tier[:raw, :cal, period])
	@bind run Select(runs, default=DataRun(0))
end

# ╔═╡ c26678c9-05e2-4be7-aa31-1fe3f0524029
md" Select Usability keys which are not used"

# ╔═╡ 31a21970-e38b-41ab-bf1f-6f92db7d4293
@bind usable MultiCheckBox([:on, :ac, :no_psd, :off], default=[:off])

# ╔═╡ dd962859-e7a8-419e-87df-3c2499f013e5
begin
	 @info "Energy calibration for period $period and run $run"

    filekeys = sort(search_disk(FileKey, l200.tier[:raw, :cal, period, run]), by = x-> x.time)
    filekey = filekeys[1]
    @info "Found filekey $filekey"
	
    chinfo = channel_info(l200, filekey) |> filterby(@pf $system == :geds && $processable)
	for u in usable
		global chinfo = chinfo |> filterby(@pf $usability != u)
	end

    sel = LegendDataManagement.ValiditySelection(filekey.time, :cal)
end;

# ╔═╡ 332edf26-da4a-4b09-a984-650c342b06e8
begin
	log_folder = joinpath(l200.tier[:log, :cal, period, run])
	log_filename = joinpath(log_folder, format("{}-{}-{}-{}-energy_calibration.md", string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
	if isfile(log_filename)
		@info "Load processing log"
		Markdown.parse_file(log_filename)
	end;
end

# ╔═╡ 9fe1dfb6-b093-4e23-92ca-71680b60f3f4
md"Select filter type"

# ╔═╡ 8c15ffad-4a3f-49da-8278-4298c11c7d4a
begin
	e_types = [ft for ft in Symbol.(l200.metadata.dataprod.config.energy(sel).default.energy_types)]
	@bind e_type Select(e_types, default=:e_trap)
end	

# ╔═╡ 634416d9-bedf-4c4c-b08d-b5ed2dae82a0
md"Select peak for resolution plot"

# ╔═╡ 4878e314-bdea-4354-960f-8f86e385e295
@bind res_peak_plot Select(Symbol.(["Tl208a", "Bi212a", "Tl208b", "Tl208DEP", "Bi212FEP", "Tl208SEP", "Tl208FEP", "Qbb"]), default=:Qbb)

# ╔═╡ 3e553316-461c-47fe-ab83-5537f53de6a9
md"Select CTC for resolution plot"

# ╔═╡ 21f8c2f9-d88e-46af-b7d7-20eaaf715426
@bind ctc_on Select([:no_ctc, :ctc], default=:ctc)

# ╔═╡ db3beed2-4c44-465b-b877-185bf2e8860d
pars = l200.par[:cal, :energy, period, run];

# ╔═╡ 343bb91e-64cd-4061-933c-f9934ceb70c2
begin
	labels = String[]
	vlines = Int[]
	xvalues = Int[]
	notworking = Int[]
	yvalues = Float64[]
	yerrvalues = Float64[]
	pvalues = Float64[]
	perrvalues = Float64[]

	e_type_ctc = e_type
	if ctc_on == :ctc
		e_type_ctc = Symbol("$(e_type)_ctc")
	end
	for s in sort(unique(chinfo.string))
	    push!(labels, format("String:{:02d}", s))
	    push!(vlines, length(labels))
	    detectors = chinfo.detector[chinfo.string .== s]
	    for det in sort(detectors)
	        push!(labels, String(det))
	        push!(xvalues, length(labels))
	        if haskey(pars, det)
				if res_peak_plot == :Qbb
	            	val = pars[det][e_type_ctc].energy.fwhm_qbb
					val_err = pars[det][e_type_ctc].energy.fwhm_qbb_err
				else
					val = pars[det][e_type_ctc].energy[res_peak_plot].fwhm
					val_err = pars[det][e_type_ctc].energy[res_peak_plot].err.fwhm
				end
	            if val isa Number && val > 0
	                push!(yvalues, val)
	                push!(yerrvalues, val_err)
	            else
	                push!(yvalues, NaN)
	                push!(yerrvalues, NaN)
	                push!(notworking, length(labels))
	            end
	        else
	            push!(yvalues, NaN)
	            push!(yerrvalues, NaN)
	            push!(notworking, length(labels))
	        end
	        # if ch in keys(pygama_dict)
	        #     val = pygama_dict[ch]
	        #     if val isa Number && val > 0
	        #         push!(pvalues, val)
	        #         push!(perrvalues, pygama_err_dict[ch])
	        #     else
	        #         push!(pvalues, NaN)
	        #         push!(perrvalues, NaN)
	        #     end
	        # else
	        #     push!(pvalues, NaN)
	        #     push!(perrvalues, NaN)
	        # end
	    end
	end
	push!(vlines, length(labels) + 1);
	
	labelcolors = fill(:black, length(labels)+1)
	labelcolors[notworking] .= :red
	labelcolors[vlines] .= :blue
	
	function multicolor_xticks!(colors0)
	    p = Plots.current()
	    xticks, xlabels = Plots.xticks(p)[1]
	    yl = Plots.ylims(p)
	    y0 = @. zero(xticks) + yl[1] - 0.02*(yl[2] - yl[1])
	    n = length(xticks)
	    colors = colors0[1:n]
	    xticks!(xticks, vcat(fill(" ",n-1), "..............."))
	    [annotate!(xi,yi,text(li,9,ci,:right, rotation = 90)) for (xi,yi,li,ci) in zip(xticks, y0, xlabels, colors)]
	    return Plots.current()
	end
end;

# ╔═╡ f7cf5e78-69b0-45e4-9e6b-b73cb33533cb
begin
	gr(fmt = :png)
	plot(xlabel = "Detector", ylabel = "FWHM at $(string(res_peak_plot)) (keV)", size = (2000,800), label = "",
	    xticks = (eachindex(labels), labels), xrotation = 90, gridalpha = 0.5, xlims = (0, length(labels) + 2), ylims = (1,5), thickness_scaling=1.5)
	vline!(vlines,color = :black, lw = 2, label = "", left_margin=10Plots.mm, bottom_margin=15Plots.mm, 
	    legendfontsize = 12, title = format("Julia Software Stack: {} Resolution ({}-{}-{}-{})", string(res_peak_plot), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
	scatter!(xvalues, yvalues, ms = 4, label = "$e_type $ctc_on average: " * format("{:.2f}keV", mean(filter(!isnan, yvalues)[filter(!isnan, yvalues) .< 5])), color = 1, yerr = yerrvalues, legend=:bottomright)
	#scatter!(xvalues, pvalues, ms = 6, label = "Pygama: " * format("{:.2f}keV", mean(filter(!isnan, pvalues))), color = 2, yerr = perrvalues)
	multicolor_xticks!(labelcolors[1:end-1])
end

# ╔═╡ 7485563e-8f49-498b-8dba-1c375a9e5fed
selected_dets = vcat(sort(chinfo.detector), [:None]);

# ╔═╡ 0e685d5c-2cfd-4f02-90a7-846b62a6426b
md"Select detector"

# ╔═╡ 6f67faca-1065-43f6-94a2-345ac74a6a6f
@bind det Select(selected_dets, default=:None)

# ╔═╡ 0e3376e3-b56f-4529-baaf-6bcb72b054f7
begin
	ch_short = chinfo.channel[findfirst(chinfo.detector .==  det)]
	ch = format("ch{}", ch_short)
	string_number = chinfo.string[findfirst(chinfo.detector .==  det)]
end;

# ╔═╡ 06fe4ab8-c5e5-4ca9-ad72-ec5b78f6e45f
begin
	if haskey(l200.metadata.dataprod.config.energy(sel), det)
		energy_config = merge(l200.metadata.dataprod.config.energy(sel).default, l200.metadata.dataprod.config.energy(sel)[det])
		@debug "Use energy config for detector $det"
	else
		energy_config = l200.metadata.dataprod.config.energy(sel).default
		@debug "Use default config for energy"
	end

	if haskey(l200.metadata.dataprod.config.qc(sel), det)
		qc_config = merge(l200.metadata.dataprod.config.qc(sel).default, l200.metadata.dataprod.config.qc(sel)[det])
		@debug "Use qc config for detector $det"
	else
		qc_config = l200.metadata.dataprod.config.qc(sel).default
		@debug "Use default qc config"
	end

	@debug "Get $e_type CT correction factor"
	if l200.par[:cal, :energy](sel)[det][e_type].ctc.fct isa PropDicts.MissingProperty
		@error "Error in $e_type CT correction factor for channel $ch"
		fct = 0.0
		# throw(ErrorException("Error in $e_type CT correction factor"))
	else
		fct = l200.par[:cal, :energy](sel)[det][e_type].ctc.fct
	end
end;

# ╔═╡ f5140e3a-2601-489a-9098-dc78b24ec0c3
begin
	@debug "Read data for $det"
	ch_filekeys = Vector{FileKey}()
	for fk in filekeys
		if !isfile(l200.tier[:dsp, fk])
			@warn "File $(basename(l200.tier[:dsp, fk])) does not exist, skip"
			continue
		end
		if !haskey(LHDataStore(l200.tier[:dsp, fk], "r"), ch)
			@warn "Channel $ch not found in $(basename(l200.tier[:dsp, fk])), skip"
			continue
		end
		push!(ch_filekeys, fk)
	end

	if isempty(ch_filekeys)
		@error "No valid filekeys found for channel $ch ($det), skip"
		throw(LoadError("$det", 154,"No filekeys found for channel $ch ($det)"))
	end

	data_ch = fast_flatten([
		LHDataStore(
			ds -> begin
				@debug "Reading from \"$(ds.data_store.filename)\""
				ds[ch][:]
			end,
			l200.tier[:dsp, fk]
		) for fk in ch_filekeys
	])
	if length(data_ch) < 50000
		@error "Not enough data points for channel $ch ($det), skip"
		throw(ErrorException("Not enough data points for channel $ch ($det)"))
	end
end;

# ╔═╡ 4b7b64ec-9689-4a0a-acc1-8f90061e5498
md""" # Energy Calibration Investigator
"""

# ╔═╡ 9d17229a-f799-4efc-9675-718153455cbd
function qc_config_input(dsp_pars::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(par[1]): $(
				Child(par[1], Slider(par[2], default=par[3], show_value=true))
			)"""
			
			for par in dsp_pars
		]
		
		md"""
		#### QC Config
		$(inputs)
		"""
	end
end;

# ╔═╡ cad476ca-99b0-48ce-8201-c25d2949e432
@bind dsp_config_slider qc_config_input([["bl_mean_sigma", (0.1:0.1:7), qc_config.blmean.sigma], ["bl_slope_sigma", (0.1:0.1:7), qc_config.blslope.sigma], ["bl_std_sigma", (0.1:0.1:7), qc_config.blsigma.sigma], ["t0_min", (0:1:120)u"µs", 42u"µs"],  ["t0_max", (0:1:120)u"µs", 52u"µs"]])

# ╔═╡ 8fbea9ab-570e-425b-999b-9619db3361ae
begin
	result_blmean, report_blmean = get_centered_gaussian_window_cut(data_ch.blmean, qc_config.blmean.min, qc_config.blmean.max, dsp_config_slider.bl_mean_sigma,; n_bins_cut=convert(Int64, round(length(data_ch)/100)), relative_cut=qc_config.blmean.relative_cut, fixed_center=false, left=true)

	blmean_qc = result_blmean.low_cut .< data_ch.blmean .< result_blmean.high_cut
	@debug format("Baseline Mean cut surrival fraction {:.2f}%", count(blmean_qc) / length(data_ch) * 100)

	result_blsigma, report_blsigma = get_centered_gaussian_window_cut(data_ch.blsigma, qc_config.blsigma.min, qc_config.blsigma.max, dsp_config_slider.bl_std_sigma,; n_bins_cut=convert(Int64, round(length(data_ch)/100)), relative_cut=qc_config.blsigma.relative_cut, fixed_center=false, left=true)

	blsigma_qc = result_blsigma.low_cut .< data_ch.blsigma .< result_blsigma.high_cut
	@debug format("Baseline Sigma cut surrival fraction {:.2f}%", count(blsigma_qc) / length(data_ch) * 100)
	
	result_blslope, report_blslope = get_centered_gaussian_window_cut(data_ch.blslope, qc_config.blslope.min*u"ns^-1", qc_config.blslope.max*u"ns^-1", dsp_config_slider.bl_slope_sigma,; n_bins_cut=convert(Int64, round(length(data_ch)/20)), relative_cut=qc_config.blslope.relative_cut);

	blslope_qc = result_blslope.low_cut .< data_ch.blslope .< result_blslope.high_cut
    @debug format("Baseline Slope cut surrival fraction {:.2f}%", count(blslope_qc) / length(data_ch) * 100)
	
	t0_qc = dsp_config_slider.t0_min .< data_ch.t0 .< dsp_config_slider.t0_max
    @debug format("t0 cut surrival fraction {:.2f}%", count(t0_qc) / length(data_ch) * 100)

	e_cut = data_ch.e_trap .> 0 .&& .!isnan.(data_ch.e_trap) .&& .!isnan.(data_ch.e_cusp) .&& .!isnan.(data_ch.e_zac)
    @debug format("Energy cut surrival fraction {:.2f}%", count(e_cut) / length(data_ch) * 100)

	inTrace_qc = .!(data_ch.inTrace_intersect .> data_ch.t0 .+ 2 .* data_ch.drift_time .&& data_ch.inTrace_n .> 1)
	@debug format("Intrace pile-up cut surrival fraction {:.2f}%", count(inTrace_qc) / length(data_ch) * 100)

	qc = blmean_qc .&& blslope_qc .&& blsigma_qc .&& t0_qc .&& inTrace_qc .&& e_cut
	@debug format("Total QC cut surrival fraction {:.2f}%", count(qc) / length(data_ch) * 100)
end;

# ╔═╡ 186004a7-44ea-4be6-bc90-611db362bb3b
begin
	plotlyjs()
	p_blmean = plot(report_blmean, xlabel="Baseline Mean (ADC)", ylabel="Counts")
	p_blsigma = plot(report_blsigma, xlabel="Baseline Std (ADC)", ylabel="Counts")
	p_blslope = plot(report_blslope, xlabel="Baseline Slope (1/ns)", ylabel="Counts")
	p_t0 = histogram(data_ch.t0, bins=:fd, xlabel="t0", ylabel="Counts", xlims=(40, 55))
	vline!([dsp_config_slider.t0_min, dsp_config_slider.t0_max], label="Cut Window", color=:green, lw=3)
	vspan!([dsp_config_slider.t0_min, dsp_config_slider.t0_max], label="Cut Window", color=:lightgreen, alpha=0.2)
	plot([p_blmean, p_blsigma, p_blslope, p_t0]..., layout=(2,2), size=(1500, 1300), legend=:outertopright, bottom_margin=50*Plots.mm, plot_title=format("{} Julia Cal. QC Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end

# ╔═╡ 2ebed1a5-16bc-4daf-b2c4-0743146d86cc
data_ch_after_qc = data_ch[qc];

# ╔═╡ 50b016cd-859e-4509-ace9-d10e59f706e7
function fit_config_input(dsp_pars::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(par[1]): $(
				Child(par[1], Slider(par[2], default=par[3], show_value=true))
			)"""
			
			for par in dsp_pars
		]
		
		md"""
		#### Fit Config
		$(inputs)
		"""
	end
end;

# ╔═╡ df803c97-7c91-4c6a-a72c-fed229fb0956
@bind th228_names MultiSelect(energy_config.th228_names)

# ╔═╡ a91b6eed-c3d0-4ccf-af70-0125d224f192
@bind fit_config fit_config_input([["n_bins", (100:100:30000), energy_config.n_bins]])

# ╔═╡ b51a549d-c229-46b0-9fdb-818b12a5724a
@bind left_windows fit_config_input([["$(peak)_left_window_size", (0:1:50), Vector{Float64}(energy_config.left_window_sizes)[energy_config.th228_names .== peak][1]] for peak in th228_names])

# ╔═╡ b28767b2-6bee-494b-9e18-5ada7f85db07
@bind right_windows fit_config_input([["$(peak)_right_window_size", (0:1:50), Vector{Float64}(energy_config.right_window_sizes)[energy_config.th228_names .== peak][1]] for peak in th228_names])

# ╔═╡ d51f1745-8505-41ad-b627-11252b046956
begin
	th228_lines = Vector{Float64}(energy_config.th228_lines)[[peak in th228_names for peak in energy_config.th228_names]]
	# th228_names_dict  = Dict{Float64, Symbol}(th228_lines .=> Symbol.(energy_config.th228_names))
	window_sizes = Vector{Tuple{Float64, Float64}}([(l,r) for (l,r) in zip(values(left_windows), values(right_windows))])
	n_bins = fit_config.n_bins
end;

# ╔═╡ a6d0ae1e-a801-4dfe-ac0a-afcb67941aca
begin
	plotlyjs()
	quantile_perc = nothing
	if !(energy_config.quantile_perc isa Number)
		quantile_perc = parse(Float64, energy_config.quantile_perc)
	else
		quantile_perc = energy_config.quantile_perc
	end
	# quantile_perc = 0.999
	result_simple, report_simple = nothing, nothing
	if ctc_on == :ctc
		result_simple, report_simple = simple_calibration(getproperty(data_ch_after_qc, e_type) .+ fct .* data_ch_after_qc.qdrift, th228_lines, window_sizes,; n_bins=n_bins, quantile_perc=quantile_perc)
	else
		result_simple, report_simple = simple_calibration(getproperty(data_ch_after_qc, e_type), th228_lines, window_sizes,; n_bins=n_bins, quantile_perc=quantile_perc)
	end
	plot(report_simple, cal=true, size=(1200, 500), title=format("{} Simple Calibration ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end

# ╔═╡ 1a9c4b93-159e-4887-bcaa-c3aba2c3c188
begin
	plotlyjs()
	result_fit, report_fit = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names)
	peak_fit_plot = plot.(values(report_fit), titleloc=:center, titlefont=font(family="monospace",halign=:center, pointsize=20), ticks=:native, right_margin=10mm, top_margin=5mm, legend=:bottomleft; show_label=true)
	for (i, p) in enumerate(peak_fit_plot)
    	xticks!(p, convert(Int, round(xlims(p)[1], digits=0)):5:convert(Int, round(xlims(p)[2], digits=0)))
    	title!(p, string(collect(keys(result_fit))[i]))
    # if i != 1
    #     plot!(showlegend=false)
    # end
	end
	plot(
	    peak_fit_plot...,
	    framestyle=:box,
	    legend=:outerright,
	    layout=(1, 7),
	    grid=true, gridalpha=0.2, gridcolor=:black, gridlinewidth=0.5,
	    xguidefont=font(family="monospace",halign=:center, pointsize=18),
	    yguidefont=font(family="monospace",halign=:center, pointsize=18),
	    xtickfontsize=10,
	    ytickfontsize=10,
	    size=(4000, 500),
	    # margins=25mm
	)
end

# ╔═╡ 89bfab7a-9f29-49e0-a313-544125367ee8
begin
	@debug "Get calibration values"
	m_cal_simple = result_simple.c
	μ = [result_fit[p].μ for p in th228_names] ./ m_cal_simple
	μ_err = [result_fit[p].err.μ for p in th228_names] ./ m_cal_simple
	m_calib, n_calib = fit_calibration(μ, th228_lines)
	@debug format("Found calibration curve: E[keV] = {:.2f} + {:.2f}*E[ADC]", n_calib, m_calib)
end	

# ╔═╡ 900e6eb3-ca58-4399-9ee5-9f43d7ed92c3
begin
	scatter(μ, th228_lines, yerror=μ_err, ms=5, color=:black, framestyle=:box, markershape= :x, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], link=:x, label="Peak Positions", xlabel="Energy (ADC)", xlabelfontsize=10, ylabel="Energy (keV)", ylabelfontsize=10, legend=:topleft, legendfontsize=8, legendfont=font(8), legendtitlefontsize=8, legendtitlefont=font(8), xlims = (0, 21000), xticks = (0:2000:22000), margin=5mm, thickness_scaling=1.5, xformatter=:plain, size=(1400, 800))
	plot!(ylims = (0, 3000), yticks = (200:200:3000), subplot=1, xlabel="", xticks = :none, bottom_margin=-4mm)
	plot!(0:1:20000, x -> m_calib* x + n_calib, label="Best Fit: $(round(n_calib, digits=2)) + x*$(round(m_calib, digits=2)))", line_width=2, color=:red, subplot=1, xformatter=_->"")
	plot!(μ, ((m_calib .* μ .+ n_calib) .- th228_lines) ./ th228_lines .* 100 , label="Residuals", ylabel="Residuals (%)", line_width=2, color=:black, st=:scatter, ylims = (-0.1, 0.1), markershape=:x, subplot=2, legend=:topleft, top_margin=0mm, framestyle=:box)
	plot!(legend = :topleft, title=format("{} Calibration Curve ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)), subplot=1)
end

# ╔═╡ 82c13dd6-9abe-4efa-8cc4-97817db0c62c
begin
	fwhm     = ([result_fit[p].fwhm for p in th228_names] ./ m_cal_simple) .* m_calib
	fwhm_err = ([result_fit[p].err.fwhm for p in th228_names] ./ m_cal_simple) .* m_calib
	result_fwhm, report_fwhm = fit_fwhm(th228_lines, fwhm)
	@debug "Found FWHM: $(round(result_fwhm.qbb, digits=2)) +- $(round(result_fwhm.err.qbb, digits=2))keV"
end

# ╔═╡ c1c66dfc-5f14-4071-98d5-b8d1b0a615a1
begin
	scatter(th228_lines, fwhm, yerror=fwhm_err, ms=5, color=:black, framestyle=:box, markershape= :x, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], link=:x, label="Peak FWHMs", xlabel="Energy (keV)", xlabelfontsize=10, ylabel="FWHM (keV)", ylabelfontsize=10, legend=:topleft, legendfontsize=8, legendfont=font(8), legendtitlefontsize=8, legendtitlefont=font(8), xlims = (0, 3000), xticks = (convert(Int, 0):300:convert(Int, round(3000, digits=0))), margin=5mm, thickness_scaling=1.5, size=(1400, 800))
	plot!(0:0.1:3000, x -> report_fwhm.f_fit(x), label="Best Fit: Sqrt($(round(report_fwhm.v[1], digits=2)) + x*$(round(report_fwhm.v[2]*100, digits=2))e-3)", line_width=2, color=:red, subplot=1, xlabel="", xticks=:none, bottom_margin=-4mm)
	hline!([result_fwhm.qbb], label="Qbb/keV: $(round(result_fwhm.qbb, digits=2))+-$(round(result_fwhm.err.qbb, digits=2))", color=:green)
	hspan!([result_fwhm.qbb - result_fwhm.err.qbb, result_fwhm.qbb + result_fwhm.err.qbb], color=:green, alpha=0.2, label="")
	plot!(th228_lines, ((report_fwhm.f_fit.(th228_lines) .- fwhm) ./ fwhm) .* 100 , label="Residuals", ylabel="Residuals (%)", line_width=2, color=:black, st=:scatter, ylims = (-10, 10), markershape=:x, legend=:topleft, subplot=2, framestyle=:box, top_margin=0mm)
	plot!(legend = :topleft, title=format("{} FWHM ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)), subplot=1)
end

# ╔═╡ Cell order:
# ╟─ecc1218c-a408-4de6-8cb5-40d29613b9e7
# ╟─98a808df-40ff-4eda-89ee-772147f7e42f
# ╟─f2cf0e3f-d544-4a0b-bcdf-d02bf9beb9d6
# ╟─2f337776-68cd-4b38-9047-88742bfa1c8a
# ╟─4933bfe7-2d34-4283-bd4f-2cb0006c5158
# ╟─d35a8320-ae1e-485d-8751-1a2b36f2b809
# ╟─24b387ff-5df5-45f9-8857-830bc6580857
# ╟─493bdebf-7802-4938-8693-5e0648bd8a2b
# ╟─39b468d4-673c-4834-b042-d80cd9e0dfe7
# ╟─2c4dd859-0aa6-4cd2-94b3-7a171368065b
# ╟─63b6be12-3b55-450a-8a2a-3d410f0dcb6c
# ╟─36232152-9811-4520-ba65-0f855e4ebbe4
# ╟─c26678c9-05e2-4be7-aa31-1fe3f0524029
# ╟─31a21970-e38b-41ab-bf1f-6f92db7d4293
# ╟─dd962859-e7a8-419e-87df-3c2499f013e5
# ╟─332edf26-da4a-4b09-a984-650c342b06e8
# ╟─9fe1dfb6-b093-4e23-92ca-71680b60f3f4
# ╟─8c15ffad-4a3f-49da-8278-4298c11c7d4a
# ╟─634416d9-bedf-4c4c-b08d-b5ed2dae82a0
# ╟─4878e314-bdea-4354-960f-8f86e385e295
# ╟─3e553316-461c-47fe-ab83-5537f53de6a9
# ╟─21f8c2f9-d88e-46af-b7d7-20eaaf715426
# ╟─db3beed2-4c44-465b-b877-185bf2e8860d
# ╟─343bb91e-64cd-4061-933c-f9934ceb70c2
# ╟─f7cf5e78-69b0-45e4-9e6b-b73cb33533cb
# ╟─7485563e-8f49-498b-8dba-1c375a9e5fed
# ╟─0e685d5c-2cfd-4f02-90a7-846b62a6426b
# ╟─6f67faca-1065-43f6-94a2-345ac74a6a6f
# ╟─0e3376e3-b56f-4529-baaf-6bcb72b054f7
# ╟─06fe4ab8-c5e5-4ca9-ad72-ec5b78f6e45f
# ╟─f5140e3a-2601-489a-9098-dc78b24ec0c3
# ╟─4b7b64ec-9689-4a0a-acc1-8f90061e5498
# ╟─cf295ad2-13ca-4508-bf51-5cfc52229fbd
# ╟─9d17229a-f799-4efc-9675-718153455cbd
# ╟─cad476ca-99b0-48ce-8201-c25d2949e432
# ╟─8fbea9ab-570e-425b-999b-9619db3361ae
# ╟─186004a7-44ea-4be6-bc90-611db362bb3b
# ╟─2ebed1a5-16bc-4daf-b2c4-0743146d86cc
# ╟─50b016cd-859e-4509-ace9-d10e59f706e7
# ╟─df803c97-7c91-4c6a-a72c-fed229fb0956
# ╟─a91b6eed-c3d0-4ccf-af70-0125d224f192
# ╟─b51a549d-c229-46b0-9fdb-818b12a5724a
# ╟─b28767b2-6bee-494b-9e18-5ada7f85db07
# ╟─d51f1745-8505-41ad-b627-11252b046956
# ╟─a6d0ae1e-a801-4dfe-ac0a-afcb67941aca
# ╟─1a9c4b93-159e-4887-bcaa-c3aba2c3c188
# ╟─89bfab7a-9f29-49e0-a313-544125367ee8
# ╟─900e6eb3-ca58-4399-9ee5-9f43d7ed92c3
# ╟─82c13dd6-9abe-4efa-8cc4-97817db0c62c
# ╟─c1c66dfc-5f14-4071-98d5-b8d1b0a615a1
