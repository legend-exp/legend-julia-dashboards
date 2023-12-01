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
	# import Pkg; Pkg.instantiate(); Pkg.precompile() # load packages
	
	using LegendDataManagement, PropertyFunctions, TypedTables, PropDicts
	using Unitful, Formatting, LaTeXStrings
	using LegendHDF5IO, LegendDSP, LegendSpecFits
	using Distributed, ProgressMeter

	using LegendDataTypes
	using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked

	
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
## Julia Analysis Software - Inspector
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
	@info "Investigate AoE filter optimization for period $period and run $run"
	
	filekeys = sort(search_disk(FileKey, l200.tier[:raw, :cal, period, run]), by = x-> x.time)
	filekey = filekeys[1]
	@info "Found filekey $filekey"
	chinfo = channel_info(l200, filekey) |> filterby(@pf $system == :geds && $processable && $usability == :on)
	
	sel = LegendDataManagement.ValiditySelection(filekey.time, :cal)
	dsp_meta = l200.metadata.dataprod.config.cal.dsp(sel).default
	dsp_config = create_dsp_config(dsp_meta)
	@debug "Loaded DSP config: $(dsp_config)"
	
	pars_tau            = l200.par[:cal, :decay_time, period, run]
	@debug "Loaded decay times"
	
	pars_optimization       = l200.par[:cal, :optimization, period, run]
	@debug "Loaded optimization parameters"
end

# ╔═╡ f776fef4-15d9-4344-a159-65ec4a3d1188
begin
	log_folder = joinpath(l200.tier[:log, :cal, period, run])
	log_filename = joinpath(log_folder, format("{}-{}-{}-{}-sg_filter_optimization.md", string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
	if isfile(log_filename)
		@info "Load processing log"
		Markdown.parse_file(log_filename)
	end;
end

# ╔═╡ 0e685d5c-2cfd-4f02-90a7-846b62a6426b
md"Select detector"

# ╔═╡ cba6bf52-a465-4905-93db-f056f7d9f9c5
selected_dets = vcat(chinfo.detector, [:None]);

# ╔═╡ 6f67faca-1065-43f6-94a2-345ac74a6a6f
@bind det confirm(Select(sort(selected_dets), default=:None))

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
	τ = 400u"µs"
	try
		τ = pars_tau[det].tau.val*u"µs"
	catch e
		@warn "No decay time available for $det"
	end
	trap_rt = 10u"µs"
	trap_ft = 4u"µs"
	try
		pars_filter = pars_optimization[det]
		# get optimal filter parameters
    	trap_rt = pars_filter.trap.rt.val*u"µs"
    	trap_ft = pars_filter.trap.ft.val*u"µs"
	catch e
		@warn "No filter optimization parameter available for $det"
		trap_rt = 10u"µs"
		trap_ft = 4u"µs"
	end
	if haskey(l200.metadata.dataprod.config.cal.dsp(sel).optimization, det)
        optimization_config = merge(l200.metadata.dataprod.config.cal.dsp(sel).optimization.default, l200.metadata.dataprod.config.cal.dsp(sel).optimization[det])
        @debug "Use config for detector $det"
    else
        optimization_config = l200.metadata.dataprod.config.cal.dsp(sel).optimization.default
        @debug "Use default config"
    end
end;

# ╔═╡ 7b56d020-5342-4a31-9902-57b77fe2c9ad
begin
	filename = joinpath(l200.tier[DataTier(:peaks), :cal, filekey.period, filekey.run], format("{}-{}-{}-{}-{}-tier_peaks.lh5", string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category), ch))
	if !isfile(filename)
		@warn "File $filename does not exist, Skip channel $ch"
		# throw(LoadError(string(basename(filename)), 154,"File $(basename(filename)) does not exist"))
	end
end;

# ╔═╡ 73d945de-e5c7-427e-a75a-b0df28f86bd4
begin
	using StatsBase, Dates
	data_ch = LHDataStore(filename, "r")

	@debug "Loading Tl208 FEP data from $(filename)"
	wvfs_ch_dep_bi121fep = data_ch[ch].Tl208DEP_Bi212FEP.waveform[:]
	ts_ch_dep_bi121fep   = unix2datetime.(ustrip.(data_ch[ch].Tl208DEP_Bi212FEP.timestamp[:]))
	e_ch_dep_bi121fep    = data_ch[ch].Tl208DEP_Bi212FEP.daqenergy[:]
	wvfs_ch_sep   	     = data_ch[ch].Tl208SEP.waveform[:]
	ts_ch_sep 			 = unix2datetime.(ustrip.(data_ch[ch].Tl208SEP.timestamp[:]))
	
	close(data_ch)
end

# ╔═╡ 41c02aee-5f40-4af7-8b22-c23efd4da269
md"Select DEP and Bi-212 FEP separation quantile"

# ╔═╡ d668effc-aecb-446e-9f40-643221041549
@bind dep_sep_quantile Slider(0:0.01:1.0, default=optimization_config.sg.dep_sep_quantile, show_value=true)

# ╔═╡ a6a9295a-e642-4e58-9552-bda1261d2027
begin
	wvfs_ch_dep   = wvfs_ch_dep_bi121fep[e_ch_dep_bi121fep .< quantile(e_ch_dep_bi121fep, dep_sep_quantile)]
	ts_ch_dep     = ts_ch_dep_bi121fep[e_ch_dep_bi121fep .< quantile(e_ch_dep_bi121fep, dep_sep_quantile)]
end;	

# ╔═╡ 4b7b64ec-9689-4a0a-acc1-8f90061e5498
md""" # Filter Optimization Investigator
"""

# ╔═╡ cdebe8d9-71e0-43cc-b5c2-3da10213282a
function qc_config_input(dsp_pars::Vector)
	
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

# ╔═╡ a1fd9d1f-57c3-49dc-8955-fdcf3a4ed3bb
@bind qc_config_slider qc_config_input([
	["dep_nbins_blslope_cut", (0:50:5000), optimization_config.sg.cuts.dep.nbins_blslope_cut],
	["sep_nbins_blslope_cut", (0:50:5000), optimization_config.sg.cuts.sep.nbins_blslope_cut],
	["dep_rel_cut_blslope_cut", (0:0.1:1.0), optimization_config.sg.cuts.dep.rel_cut_blslope_cut],
	["sep_rel_cut_blslope_cut", (0:0.1:1.0), optimization_config.sg.cuts.sep.rel_cut_blslope_cut],
	["t0_min", (0:1:120)u"µs", 42u"µs"],  ["t0_max", (0:1:120)u"µs", 52u"µs"]])

# ╔═╡ 6c9eb359-784a-48fb-9134-77ce29f6e255
begin
	stephist(e_ch_dep_bi121fep, bins=minimum(e_ch_dep_bi121fep):5:maximum(e_ch_dep_bi121fep), xlabel="Energy (ADC)", ylabel="Counts", label="FC energy", plot_title=format("{} Julia AoE Optimization Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)), size=(1500, 800), legend=:outertopright)
	vline!([quantile(e_ch_dep_bi121fep, dep_sep_quantile)], label="DEP Bi-212 separator")
end

# ╔═╡ 71f6a835-ebac-49e6-aac1-bbcc4d73bfe8
begin
	wvfs_plot_colors  = [p for p in palette(:twelvebitrainbow, 10)]
	vline_plot_colors = [p for p in palette(:tab10, 6)]
end;

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
		#### DSP Config
		$(inputs)
		"""
	end
end;

# ╔═╡ ee95f35f-4c9c-4323-8f97-4beafab379fe
@bind dsp_config_slider dsp_config_input([["bl_mean_min", (0:1:20)u"µs", dsp_config.bl_mean[1]], ["bl_mean_max", (0:1:60)u"µs", dsp_config.bl_mean[2]], ["t0_threshold", (1:1:10), dsp_config.t0_threshold], 
	["a_grid_wl_sg_step", (0:16:160)u"ns", step(dsp_config.a_grid_wl_sg)], ["a_grid_wl_sg_start", (0:2:500)u"ns", minimum(dsp_config.a_grid_wl_sg)], ["a_grid_wl_sg_stop", (10:2:500)u"ns", maximum(dsp_config.a_grid_wl_sg)]])

# ╔═╡ 03ffcd7c-9401-4237-9e1b-a3dbcc985f49
begin
	@debug "Generating DSP AoE grid for DEP data"
	using RadiationDetectorDSP, IntervalSets
	# dsp_dep = dsp_sg_optimization(wvfs_ch_dep, dsp_config, pars_tau[det].tau.val*u"µs", pars_optimization[det].trap)
	# dsp_sep = dsp_sg_optimization(wvfs_ch_sep, dsp_config, pars_tau[det].tau.val*u"µs", pars_optimization[det].trap)
	# get config parameters
	bl_mean_min, bl_mean_max    = dsp_config_slider.bl_mean_min, dsp_config_slider.bl_mean_max
	t0_threshold                = dsp_config_slider.t0_threshold
	a_grid_wl_sg                = (dsp_config_slider.a_grid_wl_sg_start:dsp_config_slider.a_grid_wl_sg_step:dsp_config_slider.a_grid_wl_sg_stop)

	# get optimal filter parameters
    rt = pars_optimization[det].trap.rt.val*u"µs"
    ft = pars_optimization[det].trap.ft.val*u"µs"

    # get baseline mean, std and slope
    bl_stats_dep = signalstats.(wvfs_ch_dep, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs_dep_bl = shift_waveform.(wvfs_ch_dep, -bl_stats_dep.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_dep_pz = deconv_flt.(wvfs_dep_bl)
    
    # t0 determination
	wvfs_dep_max = maximum.(wvfs_dep_bl.signal)
    t50_dep = LegendDSP.get_t50(wvfs_dep_pz, wvfs_dep_max)

    # get energy for filter parameters
    uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
    wvfs_flt_rtft_dep  = uflt_rtft.(wvfs_dep_pz)
    e_dep = SignalEstimator(PolynomialDNI(4, 100u"ns")).(wvfs_flt_rtft_dep, t50_dep .+ (rt + ft/2))

    # extract current with filter length in grid with second order polynominal and first derivative
    aoe_grid_dep   = ones(Float64, length(a_grid_wl_sg), length(wvfs_dep_pz))
    for (w, wl) in enumerate(a_grid_wl_sg)
        sgflt_deriv = SavitzkyGolayFilter(wl, 2, 1)
        wvfs_sgflt_deriv = sgflt_deriv.(wvfs_dep_pz)
        current_max = LegendDSP.get_wvf_maximum.(wvfs_sgflt_deriv, 20u"µs", 100u"µs")

        aoe_grid_dep[w, :]     = ustrip.(current_max) ./ e_dep
    end
	# Load DEP data and prepare Pile-up cut
    blslope_dep, t50_dep = bl_stats_dep.slope[isfinite.(e_dep)], t50_dep[isfinite.(e_dep)]
	e_dep_before_finite = e_dep
    aoe_dep, e_dep = aoe_grid_dep[:, isfinite.(e_dep)], e_dep[isfinite.(e_dep)]
end;

# ╔═╡ 899c72e9-fe30-42a7-86fd-c9cbf0fd8af6
begin
	@debug "Generating DEP QC cuts"
    # get half truncated centered cut on blslope for pile-up rejection
    result_dep_slope_cut, report_dep_slope_cut = get_centered_gaussian_window_cut(blslope_dep, -0.1u"ns^-1", 0.1u"ns^-1", 3, ; n_bins_cut=qc_config_slider.dep_nbins_blslope_cut, relative_cut=qc_config_slider.dep_rel_cut_blslope_cut)
    # Cut on blslope, energy and t0 for simple QC
    qc_cut_dep = blslope_dep .> result_dep_slope_cut.low_cut .&& blslope_dep .< result_dep_slope_cut.high_cut .&& e_dep .> optimization_config.sg.cuts.dep.min_e .&& quantile(e_dep, first(optimization_config.sg.cuts.dep.e_quantile)) .< e_dep .< quantile(e_dep, last(optimization_config.sg.cuts.dep.e_quantile)) .&& first(optimization_config.sg.cuts.dep.t50)u"µs" .< t50_dep .< last(optimization_config.sg.cuts.dep.t50)u"µs"
    aoe_dep_afterQC, e_dep_afterQC = aoe_dep[:, qc_cut_dep], e_dep[qc_cut_dep]
	@debug "DEP QC Surrival Fraction: $(round(count(qc_cut_dep)/length(qc_cut_dep)*100, digits=2))%"
end;

# ╔═╡ 036b2882-e709-431e-ac3b-6278f3179170
begin
	stephist(e_dep, bins=quantile(e_dep, first(optimization_config.sg.cuts.dep.e_quantile)):3:quantile(e_dep, last(optimization_config.sg.cuts.dep.e_quantile)), xlabel="Energy (ADC)", ylabel="Counts", label="DEP before QC", plot_title=format("{} Julia AoE Optimization Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)), size=(1500, 800), legend=:outertopright)
	# stephist(e_dep, bins=5000, xlabel="Energy (ADC)", ylabel="Counts", label="DEP after QC")
	stephist!(e_dep_afterQC, bins=quantile(e_dep, first(optimization_config.sg.cuts.dep.e_quantile)):3:quantile(e_dep, last(optimization_config.sg.cuts.dep.e_quantile)), xlabel="Energy (ADC)", ylabel="Counts", label="DEP after QC")
	# stephist!(e_ch_dep_bi121fep, bins=1000)
end

# ╔═╡ 62abd3d7-c58d-472e-a3a2-f8b709ee6262
begin
	@debug "Generating DSP AoE grid for SEP data"
	# dsp_dep = dsp_sg_optimization(wvfs_ch_dep, dsp_config, pars_tau[det].tau.val*u"µs", pars_optimization[det].trap)
	# dsp_sep = dsp_sg_optimization(wvfs_ch_sep, dsp_config, pars_tau[det].tau.val*u"µs", pars_optimization[det].trap)
	# get config parameters

    # get baseline mean, std and slope
    bl_stats_sep = signalstats.(wvfs_ch_sep, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs_sep_bl = shift_waveform.(wvfs_ch_sep, -bl_stats_sep.mean)

    # deconvolute waveform
    wvfs_sep_pz = deconv_flt.(wvfs_sep_bl)
    
    # t0 determination
	wvfs_sep_max = maximum.(wvfs_sep_bl.signal)
    t50_sep = LegendDSP.get_t50(wvfs_sep_pz, wvfs_sep_max)

    # get energy for filter parameters
    wvfs_flt_rtft_sep  = uflt_rtft.(wvfs_sep_pz)
    e_sep = SignalEstimator(PolynomialDNI(4, 100u"ns")).(wvfs_flt_rtft_sep, t50_sep .+ (rt + ft/2))

    # extract current with filter length in grid with second order polynominal and first derivative
    aoe_grid_sep   = ones(Float64, length(a_grid_wl_sg), length(wvfs_sep_pz))
    for (w, wl) in enumerate(a_grid_wl_sg)
        sgflt_deriv = SavitzkyGolayFilter(wl, 2, 1)
        wvfs_sgflt_deriv = sgflt_deriv.(wvfs_sep_pz)
        current_max = LegendDSP.get_wvf_maximum.(wvfs_sgflt_deriv, 20u"µs", 100u"µs")

        aoe_grid_sep[w, :]     = ustrip.(current_max) ./ e_sep
    end
	# Load DEP data and prepare Pile-up cut
    blslope_sep, t50_sep = bl_stats_sep.slope[isfinite.(e_sep)], t50_sep[isfinite.(e_sep)]
	e_sep_before_finite = e_sep
    aoe_sep, e_sep = aoe_grid_sep[:, isfinite.(e_sep)], e_sep[isfinite.(e_sep)]
end;

# ╔═╡ 1059ee18-47a5-4b7a-80ec-01c07cf77101
begin
	@debug "Generating SEP QC cuts"
    # get half truncated centered cut on blslope for pile-up rejection
    result_sep_slope_cut, report_sep_slope_cut = get_centered_gaussian_window_cut(blslope_sep, -0.1u"ns^-1", 0.1u"ns^-1", 3, ; n_bins_cut=qc_config_slider.sep_nbins_blslope_cut, relative_cut=qc_config_slider.sep_rel_cut_blslope_cut)
    # Cut on blslope, energy and t0 for simple QC
    qc_cut_sep = blslope_sep .> result_sep_slope_cut.low_cut .&& blslope_sep .< result_sep_slope_cut.high_cut .&& e_sep .> optimization_config.sg.cuts.sep.min_e .&& quantile(e_sep, first(optimization_config.sg.cuts.sep.e_quantile)) .< e_sep .< quantile(e_sep, last(optimization_config.sg.cuts.sep.e_quantile)) .&& first(optimization_config.sg.cuts.sep.t50)u"µs" .< t50_sep .< last(optimization_config.sg.cuts.sep.t50)u"µs"
    aoe_sep_afterQC, e_sep_afterQC = aoe_sep[:, qc_cut_sep], e_sep[qc_cut_sep]
	@debug "SEP QC Surrival Fraction: $(round(count(qc_cut_sep)/length(qc_cut_sep)*100, digits=2))%"
end;

# ╔═╡ c3f8213b-9538-44bb-900b-a0f2ee47d5a9
begin
	stephist(e_sep, bins=quantile(e_sep, first(optimization_config.sg.cuts.sep.e_quantile)):3:quantile(e_sep, last(optimization_config.sg.cuts.sep.e_quantile)), xlabel="Energy (ADC)", ylabel="Counts", label="SEP before QC", plot_title=format("{} Julia AoE Optimization Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)), size=(1500, 800), legend=:outertopright)
	# stephist(e_dep, bins=5000, xlabel="Energy (ADC)", ylabel="Counts", label="DEP after QC")
	stephist!(e_sep_afterQC, bins=quantile(e_sep, first(optimization_config.sg.cuts.sep.e_quantile)):3:quantile(e_sep, last(optimization_config.sg.cuts.sep.e_quantile)), xlabel="Energy (ADC)", ylabel="Counts", label="SEP after QC")
	# stephist!(e_ch_dep_bi121fep, bins=1000)
end

# ╔═╡ 6aa6b3fa-007a-4e4b-8819-1f2992a83d8f
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

# ╔═╡ 82c13dd6-9abe-4efa-8cc4-97817db0c62c
@bind fit_config_slider fit_config_input([["min_aoe_quantile", (0.0:0.01:1.0), optimization_config.sg.min_aoe_quantile], 
	["max_aoe_quantile", (0.0:0.01:1.0), optimization_config.sg.max_aoe_quantile],
	["min_aoe_offset", (0.0:0.01:5.0), optimization_config.sg.min_aoe_offset],
	["max_aoe_offset", (0.0:0.01:5.0), optimization_config.sg.max_aoe_offset],
	["sep_window_left", (1.0:1.0:50.0), first(optimization_config.sg.sep_window)],
	["sep_window_right", (1.0:1.0:50.0), last(optimization_config.sg.sep_window)],
	["dep_window_left", (1.0:1.0:50.0), first(optimization_config.sg.dep_window)],
	["dep_window_right", (1.0:1.0:50.0), last(optimization_config.sg.dep_window)],
	["nbins_dep_cut", (10:10:1000), optimization_config.sg.nbins_dep_cut],
	["dep_rel_cut", (0.0:0.05:1.0), optimization_config.sg.dep_rel_cut]
])

# ╔═╡ dcf89608-c309-4a6a-ade5-b1b826db6747
begin
	dep = optimization_config.sg.dep
	dep_window = [fit_config_slider.dep_window_left, fit_config_slider.dep_window_right]
    sep = optimization_config.sg.sep
	sep_window = [fit_config_slider.sep_window_left, fit_config_slider.sep_window_right]
	 
    
	# create empty arrays for sf and sf_err
    sep_sfs     = ones(length(a_grid_wl_sg)) .* 100
    sep_sfs_err = zeros(length(a_grid_wl_sg))

	reports_sep_fit = NamedTuple[]
	results_sep_fit = NamedTuple[]

	psd_cuts = []

	# prepare peakhist
    result_dep, report_dep = prepare_dep_peakhist(e_dep_afterQC, dep; n_bins_cut=fit_config_slider.nbins_dep_cut, relative_cut=fit_config_slider.dep_rel_cut)

    # get calib constant from fit on DEP peak
    e_dep_calib = e_dep_afterQC .* result_dep.m_calib
    e_sep_calib = e_sep_afterQC .* result_dep.m_calib

    for (i_aoe, wl) in enumerate(a_grid_wl_sg)
		
        aoe_dep_i = aoe_dep_afterQC[i_aoe, :][isfinite.(aoe_dep_afterQC[i_aoe, :])] ./ result_dep.m_calib
        e_dep_i   = e_dep_calib[isfinite.(aoe_dep_afterQC[i_aoe, :])]

		# prepare AoE
        max_aoe_dep_i = quantile(aoe_dep_i, fit_config_slider.max_aoe_quantile) + fit_config_slider.max_aoe_offset
        min_aoe_dep_i = quantile(aoe_dep_i, fit_config_slider.min_aoe_quantile) + fit_config_slider.min_aoe_offset
		try
            psd_cut = get_psd_cut(aoe_dep_i, e_dep_i; window=dep_window, cut_search_interval=(min_aoe_dep_i, max_aoe_dep_i))

            aoe_sep_i = aoe_sep_afterQC[i_aoe, :][isfinite.(aoe_sep_afterQC[i_aoe, :])] ./ result_dep.m_calib
			e_sep_i   = e_sep_calib[isfinite.(aoe_sep_afterQC[i_aoe, :])]

            result_sep, report_sep = get_peak_surrival_fraction(aoe_sep_i, e_sep_i, sep, sep_window, psd_cut.cut; uncertainty=true, low_e_tail=false)
            sep_sfs[i_aoe]     = result_sep.sf * 100
            sep_sfs_err[i_aoe] = result_sep.err.sf * 100
			push!(reports_sep_fit, report_sep)
			push!(results_sep_fit, result_sep)
			push!(psd_cuts, psd_cut.cut)
        catch
            @warn "Couldn't process window length $wl"
			push!(reports_sep_fit, NamedTuple())
			push!(results_sep_fit, NamedTuple())
			push!(psd_cuts, 0.0)
        end
    end

    # get minimal surrival fraction and window length
    if isempty(sep_sfs[1.0 .< sep_sfs .< 100])
        @warn "No valid SEP SF found, setting to NaN"
        min_sf = NaN
        min_sf_err = NaN
        @warn "No valid window length found, setting to default"
        wl_sg_min_sf = 100u"ns"
    else
        min_sf     = minimum(sep_sfs[1.0 .< sep_sfs .< 100])
        min_sf_err = sep_sfs_err[sep_sfs .== min_sf][1]
        wl_sg_min_sf = a_grid_wl_sg[1.0 .< sep_sfs .< 100][findmin(sep_sfs[1.0 .< sep_sfs .< 100])[2]]
    end
    # generate result and report
    result = (
        wl = wl_sg_min_sf,
        sf = min_sf,
        sf_err = min_sf_err
    )
    report = (
        wl = result.wl,
        min_sf = result.sf,
        min_sf_err = result.sf_err,
        a_grid_wl_sg = collect(a_grid_wl_sg),
        sfs = sep_sfs,
        sfs_err = sep_sfs_err
    )
 end;

# ╔═╡ 86717233-5277-460d-a933-280b28a1ee8d
begin
	i_aoe_tst = 9
	wl = a_grid_wl_sg[i_aoe_tst]
	aoe_dep_i = aoe_dep_afterQC[i_aoe_tst, :][isfinite.(aoe_dep_afterQC[i_aoe_tst, :])] ./ result_dep.m_calib
	e_dep_i   = e_dep_calib[isfinite.(aoe_dep_afterQC[i_aoe_tst, :])]
	
	# prepare AoE
	max_aoe_dep_i = quantile(aoe_dep_i, fit_config_slider.max_aoe_quantile) + fit_config_slider.max_aoe_offset
	min_aoe_dep_i = quantile(aoe_dep_i, fit_config_slider.min_aoe_quantile) + fit_config_slider.min_aoe_offset

	psd_cut = get_psd_cut(aoe_dep_i, e_dep_i; window=dep_window, cut_search_interval=(min_aoe_dep_i, max_aoe_dep_i))
end

# ╔═╡ 032b5bd5-f68a-403f-b567-3a57d94d4c98
@bind a_wl_grid_selector Select(collect(a_grid_wl_sg), default=result.wl)

# ╔═╡ 3240d0db-c448-4e00-a350-a843fe87e414
begin
	sf_plot = plot(report, title="SEP Surrival Fractions")
	ylims!(3, 30)
	i_aoe = findfirst(a_grid_wl_sg .== a_wl_grid_selector)
	sep_sf_plot = plot(reports_sep_fit[i_aoe].after, title="WL $a_wl_grid_selector - SF: $(round(results_sep_fit[i_aoe].sf * 100, digits=2)) ± $(round(results_sep_fit[i_aoe].err.sf * 100, digits=2))")
	plot!(reports_sep_fit[i_aoe].before; show_label=false)
	sep_aoe_plot = stephist(aoe_sep_afterQC[i_aoe, :][isfinite.(aoe_sep_afterQC[i_aoe, :])] ./ result_dep.m_calib, bins=0.1:1e-3:0.7, label="SEP AoE", title="WL $a_wl_grid_selector, Cut: $(round(psd_cuts[i_aoe], digits=3))")
	vline!([psd_cuts[i_aoe]], label="PSD cut", color=:red, lw=2.5)
	plot([sf_plot, sep_sf_plot, sep_aoe_plot]..., layout=(1,3), size=(2000, 800), legend=:outertopright, bottom_margin=50*Plots.mm, plot_title=format("{} Julia AoE Filter Optimization Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end

# ╔═╡ f822eb10-4a73-4506-840a-f754aa99a767
begin
	wvfs_ch_sep_accepted = wvfs_ch_sep[isfinite.(e_sep_before_finite)][qc_cut_sep][aoe_sep_afterQC[i_aoe, :] ./ result_dep.m_calib .> psd_cuts[i_aoe]]
	ts_ch_sep_accepted = ts_ch_sep[isfinite.(e_sep_before_finite)][qc_cut_sep][aoe_sep_afterQC[i_aoe, :] ./ result_dep.m_calib .> psd_cuts[i_aoe]]
	wvfs_ch_dep_accepted = wvfs_ch_dep[isfinite.(e_dep_before_finite)][qc_cut_dep][aoe_dep_afterQC[i_aoe, :] ./ result_dep.m_calib .> psd_cuts[i_aoe]]
	ts_ch_dep_accepted = ts_ch_dep[isfinite.(e_dep_before_finite)][qc_cut_dep][aoe_dep_afterQC[i_aoe, :] ./ result_dep.m_calib .> psd_cuts[i_aoe]]
	wvfs_ch_sep_rejected = wvfs_ch_sep[isfinite.(e_sep_before_finite)][qc_cut_sep][aoe_sep_afterQC[i_aoe, :] ./ result_dep.m_calib .< psd_cuts[i_aoe]]
	ts_ch_sep_rejected = ts_ch_sep[isfinite.(e_sep_before_finite)][qc_cut_sep][aoe_sep_afterQC[i_aoe, :] ./ result_dep.m_calib .< psd_cuts[i_aoe]]
	wvfs_ch_dep_rejected = wvfs_ch_dep[isfinite.(e_dep_before_finite)][qc_cut_dep][aoe_dep_afterQC[i_aoe, :] ./ result_dep.m_calib .< psd_cuts[i_aoe]]
	ts_ch_dep_rejected = ts_ch_dep[isfinite.(e_dep_before_finite)][qc_cut_dep][aoe_dep_afterQC[i_aoe, :] ./ result_dep.m_calib .< psd_cuts[i_aoe]]
end;

# ╔═╡ d8f70629-7019-48e6-b7cd-da5a0100bcda
begin
	@bind set_n_sep Slider(1:10:length(wvfs_ch_sep_rejected), show_value=true)
end

# ╔═╡ c36241c8-637e-415e-b17f-b3273faf7ba8
md" Selected SEP wvfs subset: $set_n_sep"

# ╔═╡ 1aa1c99b-acfe-4b2c-8b9a-e420b4fe512c
idx_wvfs_plot_sep = set_n_sep:set_n_sep+10-1;

# ╔═╡ edd2689c-6d75-42d9-960f-acabea398e08
begin
	@bind set_n_dep Slider(1:10:length(wvfs_ch_dep_accepted), show_value=true)
end

# ╔═╡ 14972330-9f6d-402e-ab30-bfa94281f14c
md" Selected DEP wvfs subset: $set_n_dep"

# ╔═╡ 1e7bf71a-a7c4-4765-ab53-f9771729407c
idx_wvfs_plot_dep = set_n_dep:set_n_dep+10-1;

# ╔═╡ 43fda9f9-acb4-4497-ad60-b84accda207c
begin
	plotlyjs()
	result_dep_calib, report_dep_calib = prepare_dep_peakhist(e_dep_calib, dep; n_bins_cut=fit_config_slider.nbins_dep_cut, relative_cut=fit_config_slider.dep_rel_cut)
	plot(report_dep_calib,  size=(1500, 800), title="Initial DEP fit")
	plot!([dep-first(dep_window), dep-first(dep_window)], [1, 1000], label="DEP window", color=:orange, lw=2.5, ls=:dot)
	plot!([dep+first(dep_window), dep+first(dep_window)], [1, 1000], label="DEP window", color=:orange, lw=2.5, ls=:dot, showlegend=false)
	xlims!(1570, 1620)
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
@bind detector_plot_config_slider detector_plot_config_input([["xlims_left", (0:1:120), 46], ["xlims_right", (0:1:120), 53]])

# ╔═╡ 54860a44-0823-4cb5-8958-474948138a25
begin
	p_accepted = plot(u"µs", NoUnits, size=(1500, 500), legend=:outertopright, thickness_scaling=1.0, title="Accepted SEP Events")
	plot!(wvfs_ch_sep_accepted[ifelse(last(idx_wvfs_plot_sep) >= length(wvfs_ch_sep_accepted), length(wvfs_ch_sep_accepted)-10:length(wvfs_ch_sep_accepted), idx_wvfs_plot_sep)], label=permutedims(ts_ch_sep_accepted[ifelse(last(idx_wvfs_plot_sep) >= length(wvfs_ch_sep_accepted), length(wvfs_ch_sep_accepted)-10:length(wvfs_ch_sep_accepted), idx_wvfs_plot_sep)]), color=permutedims(wvfs_plot_colors))
	xlims!(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right)
	p_rejected = plot(u"µs", NoUnits, size=(1000, 700), legend=:topright, thickness_scaling=1.0, title="Rejected SEP Events")
	plot!(wvfs_ch_sep_rejected[idx_wvfs_plot_sep], label=permutedims(ts_ch_sep_rejected[idx_wvfs_plot_sep]), color=permutedims(wvfs_plot_colors))
	xlims!(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right)
	
	plot([p_accepted, p_rejected]..., layout=(1,2), size=(1500, 600), legend=:outertopright, plot_title=format("{} Julia AoE Filter Optmimization Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end

# ╔═╡ 06d1608f-8cfc-4512-b10a-f2371d95bf10
begin
	p_accepted_dep = plot(u"µs", NoUnits, size=(1500, 500), legend=:outertopright, thickness_scaling=1.0, title="Accepted DEP Events")
	plot!(wvfs_ch_dep_accepted[idx_wvfs_plot_dep], label=permutedims(ts_ch_dep_accepted[idx_wvfs_plot_dep]), color=permutedims(wvfs_plot_colors))
	xlims!(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right)
	p_rejected_dep = plot(u"µs", NoUnits, size=(1000, 700), legend=:topright, thickness_scaling=1.0, title="Rejected DEP Events")
	plot!(wvfs_ch_dep_rejected[ifelse(last(idx_wvfs_plot_dep) >= length(wvfs_ch_dep_rejected), length(wvfs_ch_dep_rejected)-10:length(wvfs_ch_dep_rejected), idx_wvfs_plot_dep)], label=permutedims(ts_ch_dep_rejected[ifelse(last(idx_wvfs_plot_dep) >= length(wvfs_ch_dep_rejected), length(wvfs_ch_dep_rejected)-10:length(wvfs_ch_dep_rejected), idx_wvfs_plot_dep)]), color=permutedims(wvfs_plot_colors))
	xlims!(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right)
	
	plot([p_accepted_dep, p_rejected_dep]..., layout=(1,2), size=(1500, 600), legend=:outertopright, plot_title=format("{} Julia AoE Filter Optmimization Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
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
# ╟─dd962859-e7a8-419e-87df-3c2499f013e5
# ╟─f776fef4-15d9-4344-a159-65ec4a3d1188
# ╟─0e685d5c-2cfd-4f02-90a7-846b62a6426b
# ╟─cba6bf52-a465-4905-93db-f056f7d9f9c5
# ╟─6f67faca-1065-43f6-94a2-345ac74a6a6f
# ╟─f5140e3a-2601-489a-9098-dc78b24ec0c3
# ╟─7b56d020-5342-4a31-9902-57b77fe2c9ad
# ╟─41c02aee-5f40-4af7-8b22-c23efd4da269
# ╟─d668effc-aecb-446e-9f40-643221041549
# ╟─73d945de-e5c7-427e-a75a-b0df28f86bd4
# ╟─a6a9295a-e642-4e58-9552-bda1261d2027
# ╟─4b7b64ec-9689-4a0a-acc1-8f90061e5498
# ╟─cdebe8d9-71e0-43cc-b5c2-3da10213282a
# ╟─a1fd9d1f-57c3-49dc-8955-fdcf3a4ed3bb
# ╟─6c9eb359-784a-48fb-9134-77ce29f6e255
# ╟─03ffcd7c-9401-4237-9e1b-a3dbcc985f49
# ╟─899c72e9-fe30-42a7-86fd-c9cbf0fd8af6
# ╟─036b2882-e709-431e-ac3b-6278f3179170
# ╟─62abd3d7-c58d-472e-a3a2-f8b709ee6262
# ╟─1059ee18-47a5-4b7a-80ec-01c07cf77101
# ╟─c3f8213b-9538-44bb-900b-a0f2ee47d5a9
# ╠═86717233-5277-460d-a933-280b28a1ee8d
# ╠═dcf89608-c309-4a6a-ade5-b1b826db6747
# ╟─43fda9f9-acb4-4497-ad60-b84accda207c
# ╟─032b5bd5-f68a-403f-b567-3a57d94d4c98
# ╟─3240d0db-c448-4e00-a350-a843fe87e414
# ╟─cf295ad2-13ca-4508-bf51-5cfc52229fbd
# ╟─f822eb10-4a73-4506-840a-f754aa99a767
# ╟─c36241c8-637e-415e-b17f-b3273faf7ba8
# ╟─d8f70629-7019-48e6-b7cd-da5a0100bcda
# ╟─1aa1c99b-acfe-4b2c-8b9a-e420b4fe512c
# ╟─54860a44-0823-4cb5-8958-474948138a25
# ╟─14972330-9f6d-402e-ab30-bfa94281f14c
# ╟─edd2689c-6d75-42d9-960f-acabea398e08
# ╟─1e7bf71a-a7c4-4765-ab53-f9771729407c
# ╟─06d1608f-8cfc-4512-b10a-f2371d95bf10
# ╟─71f6a835-ebac-49e6-aac1-bbcc4d73bfe8
# ╟─89bfab7a-9f29-49e0-a313-544125367ee8
# ╟─ee95f35f-4c9c-4323-8f97-4beafab379fe
# ╟─6aa6b3fa-007a-4e4b-8819-1f2992a83d8f
# ╟─82c13dd6-9abe-4efa-8cc4-97817db0c62c
# ╟─cc99755f-be83-4525-8850-a9df59eaca15
# ╟─d4d256f0-45c5-4f2b-bbc2-3f44b2802805
