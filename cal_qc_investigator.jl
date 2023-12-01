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
	using TypedTables
	
	using LegendDataTypes
	using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked

	using StatsBase, Dates

	using RadiationDetectorSignals
	
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
## Julia Analysis Software - Calibration QC Investigator
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

# ╔═╡ dd962859-e7a8-419e-87df-3c2499f013e5
begin
	@info "Investigate DSP for period $period and run $run"
	
	filekeys = sort(search_disk(FileKey, l200.tier[:raw, :cal, period, run]), by = x-> x.time)
	filekey = filekeys[1]
	@info "Found filekey $filekey"
	chinfo = channel_info(l200, filekey) |> filterby(@pf $system == :geds && $processable)
	
	sel = LegendDataManagement.ValiditySelection(filekey.time, :cal)
	dsp_meta = l200.metadata.dataprod.config.cal.dsp(sel).default
	dsp_config = create_dsp_config(dsp_meta)
	@debug "Loaded DSP config: $(dsp_config)"
	
	pars_tau = l200.par[:cal, :decay_time, period, run]
	@debug "Loaded decay times"
	
	pars_optimization = l200.par[:cal, :optimization, period, run] 
	@debug "Loaded optimization parameters"
end

# ╔═╡ 0e685d5c-2cfd-4f02-90a7-846b62a6426b
md"Select detector"

# ╔═╡ 6f67faca-1065-43f6-94a2-345ac74a6a6f
@bind det Select(sort(chinfo.detector), default=:V09372A)

# ╔═╡ f5140e3a-2601-489a-9098-dc78b24ec0c3
begin
	i = findfirst(x -> x == det, chinfo.detector)
	ch_short = chinfo.channel[i]
	ch = format("ch{}", ch_short)
	string_number = chinfo.string[i]
	pos_number = chinfo.position[i]
	cc4_name = chinfo.cc4[i] * "$(chinfo.cc4ch[i])"
	# check if channel can be processed
	if !haskey(pars_tau, det)
    	@warn "No decay time for detector $det, skip channel $ch"
	else
		@info "Selected detector: $det at string $string_number position $pos_number" 
	end
	τ = pars_tau[det].tau.val*u"µs"
	pars_filter = pars_optimization[det]
end;

# ╔═╡ 73d945de-e5c7-427e-a75a-b0df28f86bd4
begin
	if haskey(l200.metadata.dataprod.config.cal.qc(sel), det)
    	qc_config = merge(l200.metadata.dataprod.config.cal.qc(sel).default, l200.metadata.dataprod.config.cal.qc(sel)[det])
    	@debug "Use config for detector $det"
	else
    	qc_config = l200.metadata.dataprod.config.cal.qc(sel).default
    	@debug "Use default config"
	end
	ch_filekeys = Vector{FileKey}()
	ch_filekeys_idx = TypedTables.Table(fk = Vector{FileKey}(), fk_idx = Vector{Int64}())
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
		append!(ch_filekeys_idx.fk, fill(fk, length(LHDataStore(l200.tier[:dsp, fk], "r")[ch])))
		append!(ch_filekeys_idx.fk_idx, collect(eachindex(LHDataStore(l200.tier[:dsp, fk], "r")[ch])))
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
	
	@debug "Number of events: $(length(data_ch))"
	@debug "Median of baseline mean : $(median(data_ch.blmean))"
	@debug "Median of baseline std  : $(median(data_ch.blsigma))"
	@debug "Median of baseline slope: $(median(data_ch.blslope))"
end;	

# ╔═╡ 60f797df-0042-45c9-bb22-3d6226677bfc
begin
	@bind selected_filekeys confirm(MultiSelect(ch_filekeys))
end

# ╔═╡ 7d65085c-fa99-4272-bb2a-b1c20650fbba
data_raw_ch = [
    LHDataStore(
        ds -> begin
            @debug "Reading from \"$(ds.data_store.filename)\""
            ds[ch*"/raw/"][:]
        end,
        l200.tier[:raw, fk]
	) for fk in selected_filekeys];

# ╔═╡ 4b7b64ec-9689-4a0a-acc1-8f90061e5498
md""" # DSP investigator
"""

# ╔═╡ 89bfab7a-9f29-49e0-a313-544125367ee8
function dsp_config_input(dsp_pars::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(par[1]): $(
				Child(par[1], Slider(par[2], default=par[3], show_value=true))
			)"""
			
			for par in dsp_pars
		]
		
		md"""
		#### Cut Config
		$(inputs)
		"""
	end
end;

# ╔═╡ ee95f35f-4c9c-4323-8f97-4beafab379fe
@bind dsp_config_slider dsp_config_input([["bl_mean_sigma", (0.1:0.1:7), qc_config.blmean.sigma], ["bl_slope_sigma", (0.1:0.1:7), qc_config.blslope.sigma], ["bl_std_sigma", (0.1:0.1:7), qc_config.blsigma.sigma], ["t0_min", (0:1:120)u"µs", 42u"µs"],  ["t0_max", (0:1:120)u"µs", 52u"µs"]])

# ╔═╡ 35fc04da-4adb-4af6-8f01-452b68577f87
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

	e_cut = data_ch.e_trap .> 0 .&& .!isnan.(data_ch.e_trap)
    @debug format("Energy cut surrival fraction {:.2f}%", count(e_cut) / length(data_ch) * 100)

	inTrace_qc = .!(data_ch.inTrace_intersect .> data_ch.t0 .+ 2 .* data_ch.drift_time .&& data_ch.inTrace_n .> 1)
	@debug format("Intrace pile-up cut surrival fraction {:.2f}%", count(inTrace_qc) / length(data_ch) * 100)

	qc = blmean_qc .&& blslope_qc .&& blsigma_qc .&& t0_qc .&& inTrace_qc .&& e_cut
	@debug format("Total QC cut surrival fraction {:.2f}%", count(qc) / length(data_ch) * 100)
end;

# ╔═╡ e2e070a1-0517-43d3-b94e-f38380a7012e
begin
	p_blmean = plot(report_blmean, xlabel="Baseline Mean (ADC)", ylabel="Counts")
	p_blsigma = plot(report_blsigma, xlabel="Baseline Std (ADC)", ylabel="Counts")
	p_blslope = plot(report_blslope, xlabel="Baseline Slope (1/ns)", ylabel="Counts")
	p_t0 = histogram(data_ch.t0, bins=:fd, xlabel="t0", ylabel="Counts", xlims=(40, 55))
	vline!([dsp_config_slider.t0_min, dsp_config_slider.t0_max], label="Cut Window", color=:green, lw=3)
	vspan!([dsp_config_slider.t0_min, dsp_config_slider.t0_max], label="Cut Window", color=:lightgreen, alpha=0.2)
	plot([p_blmean, p_blsigma, p_blslope, p_t0]..., layout=(2,2), size=(1500, 1300), legend=:outertopright, bottom_margin=50*Plots.mm, plot_title=format("{} Julia Cal. QC Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end

# ╔═╡ cc99755f-be83-4525-8850-a9df59eaca15
function detector_plot_config_input(dsp_pars::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(par[1]): $(
				Child(par[1], Slider(par[2], default=par[3], show_value=true))
			)"""
			
			for par in dsp_pars
		]
		
		md"""
		#### Detector Plot Config
		$(inputs)
		"""
	end
end;

# ╔═╡ d4d256f0-45c5-4f2b-bbc2-3f44b2802805
@bind detector_plot_config_slider detector_plot_config_input([["n_bins_e", (1000:200:15000), 10000], ["e_cut_left", (0:1:maximum(data_ch.e_trap[e_cut])), 0], ["e_cut_right", (0:1:maximum(data_ch.e_trap[e_cut])), 1.01*quantile(data_ch.e_trap[e_cut], 0.99)]])

# ╔═╡ 33face90-16fa-4865-918f-c7a547eec7a3
begin
	stephist(data_ch.e_trap[e_cut], bins=detector_plot_config_slider.n_bins_e, xlabel="Energy", ylabel="Counts", yscale=:log10, size=(1200, 600), label="e_trap")
	plot!(fill(detector_plot_config_slider.e_cut_left, 2), [0.1, 1e4], label="Energy Cut Window", color=:red, lw=2.5, ls=:dot)
	plot!(fill(detector_plot_config_slider.e_cut_right, 2), [0.1, 1e4], label="Energy Cut Window", color=:red, lw=2.5, ls=:dot, showlegend=false)
	xlims!(-300, 1.3*quantile(data_ch.e_trap[e_cut], 0.99))
end

# ╔═╡ 73cbfa7e-b18c-4e2f-a39f-ae68bbf4a6fb
cut_options = Dict(["BlMean", "BlStd", "BlSlope", "t0", "pile-up"] .=> [blmean_qc, blsigma_qc, blslope_qc, t0_qc, inTrace_qc]);

# ╔═╡ e6a3a08d-9c64-451f-931c-54d8c3e3d6c5
@bind cut_selector MultiCheckBox(collect(keys(cut_options)), default=["BlMean", "BlStd", "BlSlope", "t0", "pile-up"])

# ╔═╡ 56faaa8d-50e6-4b8e-aa7b-a8afd38e322d
begin
	n_load = 1000
	accepted_idx = e_cut
	for cut in cut_selector
		accepted_idx = accepted_idx .&& cut_options[cut]
	end
	fk_accepted = ch_filekeys_idx[accepted_idx .&& detector_plot_config_slider.e_cut_left .< data_ch.e_trap .< detector_plot_config_slider.e_cut_right]
	fk_rejected = ch_filekeys_idx[.!accepted_idx .&& detector_plot_config_slider.e_cut_left .< data_ch.e_trap .< detector_plot_config_slider.e_cut_right]
	@info format("Surrival fraction in energy window: {:.2f}%", length(fk_accepted)/(length(fk_accepted) + length(fk_rejected))*100)
end;

# ╔═╡ 42bc464c-d84f-42f2-864a-9904c02bb9d0
begin
	data_rejected = fast_flatten([data_raw_ch[i][filter(row -> row.fk == fk_raw, fk_rejected).fk_idx] for (i, fk_raw) in enumerate(selected_filekeys)])
	wvfs_rejected, ts_rejected = data_rejected.waveform, unix2datetime.(ustrip.(data_rejected.timestamp))
	data_accepted = fast_flatten([data_raw_ch[i][filter(row -> row.fk == fk_raw, fk_accepted).fk_idx] for (i, fk_raw) in enumerate(selected_filekeys)])
	wvfs_accepted, ts_accepted = data_accepted.waveform, unix2datetime.(ustrip.(data_accepted.timestamp))
end;

# ╔═╡ 7eec6fe1-4cb9-4efd-91a1-a905a61fa874
@bind set_n Slider(1:10:minimum([length(data_rejected), length(data_accepted)]), show_value=true)

# ╔═╡ 77ddc928-4b62-4241-8001-01ff84969272
md" Selected subset: $set_n"

# ╔═╡ 4c415d15-5dff-46bb-ba11-d1f7224546de
idx_wvfs_plot = set_n:set_n+10-1;

# ╔═╡ 71f6a835-ebac-49e6-aac1-bbcc4d73bfe8
begin
	wvfs_plot_colors  = [p for p in palette(:twelvebitrainbow, 10)]
	vline_plot_colors = [p for p in palette(:tab10, 6)]
end;

# ╔═╡ 97f7da06-a687-4c62-a8d3-bc42430a9ff1
begin
	p_accepted = plot(u"µs", NoUnits, size=(1500, 500), legend=:outertopright, thickness_scaling=1.0, title="Accepted Events")
	plot!(wvfs_accepted[idx_wvfs_plot], label=permutedims(ts_accepted[idx_wvfs_plot]), color=permutedims(wvfs_plot_colors))
	# p_rejected = plot(u"µs", NoUnits, size=(1000, 700), legend=:topright, thickness_scaling=1.0, title="Rejected Events")
	# plot!(wvfs_rejected[idx_wvfs_plot], label=permutedims(ts_rejected[idx_wvfs_plot]), color=permutedims(wvfs_plot_colors))
	
	# plot([p_accepted, p_rejected]..., layout=(1,2), size=(1500, 500), legend=:outertopright, plot_title=format("{} Julia Cal QC Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end

# ╔═╡ 64f2fef7-57f7-4c7a-ba83-99be01e3082b
begin
	p_rejected = plot(u"µs", NoUnits, size=(1500, 500), legend=:outertopright, thickness_scaling=1.0, title="Rejected Events")
	plot!(wvfs_rejected[idx_wvfs_plot], label=permutedims(ts_rejected[idx_wvfs_plot]), color=permutedims(wvfs_plot_colors))
end

# ╔═╡ 43204bc4-c869-4a76-ba93-500a33c39b89
# vline_options = Dict(["t0", "t10", "t50", "t90", "t99", "area1", "area2"] .=> [t0, t10, t50, t90, t99, t0 .+ dsp_config_slider.qdrift_window, t0 .+ (2 * dsp_config_slider.qdrift_window)]);

# ╔═╡ b5ddd0ba-a771-4f2e-94e1-ab9e994d69b3
# @bind vline_type_selector MultiCheckBox(sort(collect(keys(vline_options))))

# ╔═╡ bb81cf12-32c9-4bd0-92a8-b727fcf9b098
function detector_config_input(dsp_pars::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(par[1]): $(
				Child(par[1], Slider(par[2], default=par[3], show_value=true))
			)"""
			
			for par in dsp_pars
		]
		
		md"""
		#### Detector config
		$(inputs)
		"""
	end
end;

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
# ╟─dd962859-e7a8-419e-87df-3c2499f013e5
# ╟─0e685d5c-2cfd-4f02-90a7-846b62a6426b
# ╟─6f67faca-1065-43f6-94a2-345ac74a6a6f
# ╟─f5140e3a-2601-489a-9098-dc78b24ec0c3
# ╟─73d945de-e5c7-427e-a75a-b0df28f86bd4
# ╟─60f797df-0042-45c9-bb22-3d6226677bfc
# ╟─7d65085c-fa99-4272-bb2a-b1c20650fbba
# ╟─4b7b64ec-9689-4a0a-acc1-8f90061e5498
# ╟─cf295ad2-13ca-4508-bf51-5cfc52229fbd
# ╟─35fc04da-4adb-4af6-8f01-452b68577f87
# ╟─89bfab7a-9f29-49e0-a313-544125367ee8
# ╟─ee95f35f-4c9c-4323-8f97-4beafab379fe
# ╟─e2e070a1-0517-43d3-b94e-f38380a7012e
# ╟─cc99755f-be83-4525-8850-a9df59eaca15
# ╟─d4d256f0-45c5-4f2b-bbc2-3f44b2802805
# ╟─33face90-16fa-4865-918f-c7a547eec7a3
# ╟─73cbfa7e-b18c-4e2f-a39f-ae68bbf4a6fb
# ╟─e6a3a08d-9c64-451f-931c-54d8c3e3d6c5
# ╟─56faaa8d-50e6-4b8e-aa7b-a8afd38e322d
# ╟─42bc464c-d84f-42f2-864a-9904c02bb9d0
# ╟─77ddc928-4b62-4241-8001-01ff84969272
# ╟─7eec6fe1-4cb9-4efd-91a1-a905a61fa874
# ╟─4c415d15-5dff-46bb-ba11-d1f7224546de
# ╟─97f7da06-a687-4c62-a8d3-bc42430a9ff1
# ╟─64f2fef7-57f7-4c7a-ba83-99be01e3082b
# ╟─71f6a835-ebac-49e6-aac1-bbcc4d73bfe8
# ╟─43204bc4-c869-4a76-ba93-500a33c39b89
# ╟─b5ddd0ba-a771-4f2e-94e1-ab9e994d69b3
# ╟─bb81cf12-32c9-4bd0-92a8-b727fcf9b098
