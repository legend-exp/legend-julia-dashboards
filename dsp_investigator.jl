### A Pluto.jl notebook ###
# v0.19.38

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
import Pkg; Pkg.activate("/home/iwsatlas1/henkes/l200/cm2023/")

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
	Pkg.instantiate(); Pkg.precompile() # load packages
	
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
	using PlutoPlotly
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
## Julia Analysis Software - DSP Inspector
This tool lets you investigate a certain set of waveforms for a specififc detector and how certain parameters impact the DSP.
"""

# ╔═╡ 24b387ff-5df5-45f9-8857-830bc6580857
md""" ### Parameter settings
"""

# ╔═╡ a11bc8f3-c095-4100-96af-1e9d3eed0215
md"Select category"

# ╔═╡ 8c632350-3fb5-462a-95f9-9b26306655fa
@bind category Select([:cal, :phy], default=:cal)

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
	
	filekeys = sort(search_disk(FileKey, l200.tier[:raw, category, period, run]), by = x-> x.time)
	filekey = filekeys[1]
	@info "Found filekey $filekey"
	chinfo = channel_info(l200, filekey) |> filterby(@pf $system == :geds && $processable && $usability != :off)
	
	sel = LegendDataManagement.ValiditySelection(filekey.time, category)
	dsp_meta = l200.metadata.dataprod.config.dsp(sel).default
	dsp_config = create_dsp_config(dsp_meta)
	@debug "Loaded DSP config: $(dsp_config)"
	
	# pars_tau_folder     = joinpath(l200.tier[:par, :cal, period, run], "decay_time")
	# pars_filename       = format("{}-{}-{}-{}-decay_time.json", string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category))
	pars_tau            = l200.par[:cal, :decay_time, period, run]
	@debug "Loaded decay times"
	
	# pars_optimization_folder = joinpath(l200.tier[:par, :cal, period, run], "optimization")
	# pars_filename           = format("{}-{}-{}-{}-filter_optimization.json", string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category))
	pars_optimization       = l200.par[:cal, :optimization, period, run]
	@debug "Loaded optimization parameters"
end

# ╔═╡ 9edf7ca3-db5f-422d-ac11-f0e6fd17c087
begin
	log_folder = joinpath(l200.tier[:log, category, period, run])
	log_filename = joinpath(log_folder, format("{}-{}-{}-{}-dsp.md", string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
	if isfile(log_filename)
		@info "Load processing log"
		Markdown.parse_file(log_filename)
	end;
end

# ╔═╡ 0e685d5c-2cfd-4f02-90a7-846b62a6426b
md"Select detector"

# ╔═╡ 5df81bc3-f86f-4a85-927c-ef9f3285dc97
selected_dets = vcat(sort(chinfo.detector), [:None]);

# ╔═╡ 6f67faca-1065-43f6-94a2-345ac74a6a6f
@bind det Select(selected_dets, default=:None)

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

# ╔═╡ de7861dc-5665-4e88-ac00-af15ad1ada5c
md"Select $category Filekeys (*Multi-Select possible*)"

# ╔═╡ d8fef6df-d831-4248-aeae-84869714f76d
begin
	@bind selected_filekeys confirm(MultiSelect(filekeys))
end

# ╔═╡ 73d945de-e5c7-427e-a75a-b0df28f86bd4
if !isempty(selected_filekeys)
	data_ch = fast_flatten([
    LHDataStore(
        ds -> begin
            @debug "Reading from \"$(ds.data_store.filename)\""
            ds[ch*"/raw/"][:]
        end,
        l200.tier[:raw, fk]
	) for fk in selected_filekeys ])
end;

# ╔═╡ 4b7b64ec-9689-4a0a-acc1-8f90061e5498
md""" # DSP investigator
"""

# ╔═╡ f7211685-e046-47cc-8261-b0c54b466469
@bind n_wvfs_plot confirm(Slider(1:50, default=10))

# ╔═╡ 7ff8481d-4658-4a3f-a671-858edcf389cf
md"Number of plotted waveforms: $(n_wvfs_plot)"

# ╔═╡ 71f6a835-ebac-49e6-aac1-bbcc4d73bfe8
begin
	wvfs_plot_colors  = [p for p in palette(:twelvebitrainbow, n_wvfs_plot)]
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
@bind dsp_config_slider dsp_config_input([["bl_mean_min", (0:1:20)u"µs", dsp_config.bl_mean[1]], ["bl_mean_max", (0:1:60)u"µs", dsp_config.bl_mean[2]], ["pz_fit_min", (40:1:120)u"µs", dsp_config.pz_fit[1]], ["pz_fit_max", (40:1:120)u"µs", dsp_config.pz_fit[2]], ["t0_threshold", (1:1:10), dsp_config.t0_threshold], ["inTraceCut_std_threshold", (1:1:10), dsp_config.inTraceCut_std_threshold], ["qdrift_window", (1.5:0.1:5.0)u"µs", 2.5*u"µs"]])

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

# ╔═╡ 60f190de-8ee1-4c81-a73d-6c5f1b52632e
@bind detector_config_slider detector_config_input([["trap_rt", (1:0.1:16)u"µs", pars_filter[:trap].rt.val*u"µs"], ["trap_ft", (0.1:0.1:7)u"µs", pars_filter[:trap].ft.val*u"µs"], 
	["cusp_rt", (1:0.1:16)u"µs", pars_filter[:cusp].rt.val*u"µs"], ["cusp_ft", (0.1:0.1:7)u"µs", pars_filter[:cusp].ft.val*u"µs"], 
	["flt_length_cusp", (11.0:1.0:120.0)u"µs", dsp_config.flt_length_cusp],
	["zac_rt", (1:0.1:16)u"µs", pars_filter[:zac].rt.val*u"µs"], ["zac_ft", (0.1:0.1:7)u"µs", pars_filter[:zac].ft.val*u"µs"],
	["flt_length_zac", (11.0:1.0:120.0)u"µs", dsp_config.flt_length_zac],
	["sg_wl", (10:10:300)u"ns", 180u"ns"]])

# ╔═╡ 35fc04da-4adb-4af6-8f01-452b68577f87
begin
	using RadiationDetectorDSP, ArraysOfArrays, IntervalSets, Dates
	# get config parameters
	bl_mean_min, bl_mean_max    = dsp_config_slider.bl_mean_min, dsp_config_slider.bl_mean_max
	t0_threshold                = dsp_config_slider.t0_threshold
	pz_fit_min, pz_fit_max      = dsp_config_slider.pz_fit_min, dsp_config_slider.pz_fit_max
	inTraceCut_std_threshold    = dsp_config_slider.inTraceCut_std_threshold
	
	# get optimal filter parameters
	trap_rt = detector_config_slider.trap_rt
	trap_ft = detector_config_slider.trap_ft
	cusp_rt = detector_config_slider.cusp_rt
	cusp_ft = detector_config_slider.cusp_ft
	zac_rt = detector_config_slider.zac_rt
	zac_ft = detector_config_slider.zac_ft
	sg_wl   = detector_config_slider.sg_wl
	
	# get waveform data 
	wvfs = data_ch.waveform
	blfc = data_ch.baseline
	ts   = unix2datetime.(ustrip.(data_ch.timestamp))
	evID = data_ch.eventnumber
	efc  = data_ch.daqenergy

	# get CUSP and ZAC filter length and flt scale
    flt_length_zac              = detector_config_slider.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))
    flt_length_cusp             = detector_config_slider.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))
	
	# set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"

	
	# get number of samples the waveform is saturated at low and high of FADC range
	bit_depth = 16 # of FlashCam FADC
	sat_low, sat_high = 0, 2^bit_depth - bit_depth
	sat_stats = saturation.(wvfs, sat_low, sat_high)
	
	# get baseline mean, std and slope
	bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)
	
	# pretrace difference 
	pretrace_diff = flatview(wvfs.signal)[1, :] - bl_stats.mean
	
	# substract baseline from waveforms
	wvfs_fc = wvfs
	wvfs = shift_waveform.(wvfs, -bl_stats.mean)
	
	# extract decay times
	tail_stats = tailstats.(wvfs, pz_fit_min, pz_fit_max)
	
	# deconvolute waveform
	deconv_flt = InvCRFilter(τ)
	wvfs_pz = deconv_flt.(wvfs)
	
	# get tail mean, std and slope
	pz_stats = signalstats.(wvfs_pz, pz_fit_min, pz_fit_max)
	
	# get wvf maximum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)

    # t0 determination
    t0 = LegendDSP.get_t0(wvfs_pz, t0_threshold)

    # t50 determination
    t50 = LegendDSP.get_t50(wvfs_pz, wvf_max)

    # t80 determination
    t80 = LegendDSP.get_t80(wvfs_pz, wvf_max)

	# sanity --> replace NaNs to have consistency with the other code
    replace!(t0, NaN*unit(t0[1]) => zero(t0[1]))
    replace!(t50, NaN*unit(t50[1]) => zero(t50[1]))
    replace!(t80, NaN*unit(t80[1]) => zero(t80[1]))
	
	# get risetimes and drift times by intersection
	flt_intersec_90RT = Intersect(mintot = 100u"ns")
	flt_intersec_99RT = Intersect(mintot = 20u"ns")
	flt_intersec_lowRT = Intersect(mintot = 600u"ns")
	
	wvf_max = maximum.(wvfs.signal)
	t10 = flt_intersec_lowRT.(wvfs_pz, wvf_max .* 0.1).x
	t50 = flt_intersec_lowRT.(wvfs_pz, wvf_max .* 0.5).x
	t90 = flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.9).x
	t99 = flt_intersec_99RT.(wvfs_pz, wvf_max .* 0.99).x
	
	rt1090     = uconvert.(u"ns", flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.9).x - flt_intersec_lowRT.(wvfs_pz, wvf_max .* 0.1).x)
	rt1099     = uconvert.(u"ns", flt_intersec_99RT.(wvfs_pz, wvf_max .* 0.99).x - flt_intersec_lowRT.(wvfs_pz, wvf_max .* 0.1).x)
	rt9099     = uconvert.(u"ns", flt_intersec_99RT.(wvfs_pz, wvf_max .* 0.99).x - flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.90).x)
	drift_time = uconvert.(u"ns", flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.90).x - t0)
	
	# get Q-drift parameter
	int_flt = IntegratorFilter(1)
	wvfs_flt_int = int_flt.(wvfs_pz)

	area1 = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0 .+ dsp_config_slider.qdrift_window) .- SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0)
	area2 = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0 .+ (2 * dsp_config_slider.qdrift_window))   .- SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0 .+ dsp_config_slider.qdrift_window)
	qdrift = area2 .- area1

	# get LQ parameter
    area1_lq = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t80 .+ 2.5u"µs") .- SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t80)
    area2_lq = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t80 .+ 5u"µs")   .- SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t80 .+ 2.5u"µs")
    lq    = area2 .- area1
	
	# extract energy and ENC noise param from maximum of filtered wvfs
	uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
	
	wvfs_flt_10410 = uflt_10410.(wvfs_pz)
	e_10410        = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_10410, t50 .+ 12u"µs")

	# extract energy and ENC noise param from maximum of filtered wvfs
    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")

    wvfs_flt = uflt_313.(wvfs_pz)
    e_313  = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt, t50 .+ 3.5u"µs")
	
       
    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale)

    wvfs_flt_cusp = uflt_cusp_rtft.(wvfs_pz)
    e_cusp   = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_cusp, t50 .+ (flt_length_cusp /2))

    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale)

    wvfs_flt_zac = uflt_zac_rtft.(wvfs_pz)
    e_zac    = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_zac, t50 .+ (flt_length_zac /2))
	
	
	# get energy of optimized rise and flat-top time
	uflt_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
	
	wvfs_flt_trap  = uflt_rtft.(wvfs_pz)
	e_rtft         = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_trap, t0 .+ (trap_rt + trap_ft/2))
	
	
	# extract current with filter length of 180ns with second order polynominal and first derivative
	sgflt_deriv = SavitzkyGolayFilter(sg_wl, 2, 1)
	wvfs_sgflt_deriv = sgflt_deriv.(wvfs_pz)
	current_max = maximum.(wvfs_sgflt_deriv.signal)
	
	# in-trace pile-up rejector
	flt_intersec_inTrace = Intersect(mintot = 100u"ns")
	deriv_stats = signalstats.(wvfs_sgflt_deriv, bl_mean_min + wvfs_sgflt_deriv.time[1][1], bl_mean_max)
	# threshold over std
	inTraceCut = inTraceCut_std_threshold .* deriv_stats.sigma
	
	# get position and multiplicity of in-trace pile-up
	inTrace_pileUp      = flt_intersec_inTrace.(reverse_waveform.(wvfs_sgflt_deriv), inTraceCut)
	inTrace_intersect   = wvfs_pz.time[1][end] .- inTrace_pileUp.x
	inTrace_n           = inTrace_pileUp.multiplicity
	
	# get position of current rise start
    t50_current = LegendDSP.get_t50(wvfs_sgflt_deriv, maximum.(wvfs_sgflt_deriv.signal))

    # invert waveform for DC tagging
    wvfs_pz_inv = multiply_waveform.(wvfs_pz, -1.0)

    # get inverted waveform maximum
    wvfs_flt_10410_inv = uflt_10410.(wvfs_pz_inv)
    e_10410_max_inv  = maximum.(wvfs_flt_10410_inv.signal)

    wvfs_flt_313_inv = uflt_313.(wvfs_pz_inv)
    e_313_max_inv  = maximum.(wvfs_flt_10410_inv.signal)

    # t0 determination
    t0_inv = LegendDSP.get_t0(wvfs_pz_inv, t0_threshold)
end;

# ╔═╡ 73cbfa7e-b18c-4e2f-a39f-ae68bbf4a6fb
wvfs_options = Dict(["Wvf", "Wvf Bl", "Wvf Pz", "Wvf Trap", "Wvf CUSP", "Wvf ZAC", "Wvf deriv"] .=> [wvfs_fc, wvfs, wvfs_pz, wvfs_flt_trap, wvfs_flt_cusp, wvfs_flt_zac, wvfs_sgflt_deriv]);

# ╔═╡ e6a3a08d-9c64-451f-931c-54d8c3e3d6c5
@bind wvfs_type_selector MultiCheckBox(collect(keys(wvfs_options)), default=["Wvf Bl"])

# ╔═╡ 43204bc4-c869-4a76-ba93-500a33c39b89
vline_options = Dict(["t0", "t10", "t50", "t90", "t99", "area1", "area2", "ftp_trap", "ftp_zac/cusp"] .=> [t0, t10, t50, t90, t99, t0 .+ dsp_config_slider.qdrift_window, t0 .+ (2 * dsp_config_slider.qdrift_window), t50 .+ (trap_rt + trap_ft/2), t50 .+ (flt_length_cusp /2)]);

# ╔═╡ b5ddd0ba-a771-4f2e-94e1-ab9e994d69b3
@bind vline_type_selector MultiCheckBox(sort(collect(keys(vline_options))))

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
@bind detector_plot_config_slider detector_plot_config_input([["n_bins_efc", (200:200:8000), 2000], ["bin_width_e", (0.0:0.1:20.0), 10.0], ["efc_cut_left", (0:1:maximum(efc)), 0], ["efc_cut_right", (0:1:maximum(efc)), maximum(efc)], ["xlims_left", (0:1:120), 0], ["xlims_right", (0:1:120), 120]])

# ╔═╡ 56faaa8d-50e6-4b8e-aa7b-a8afd38e322d
begin
	using Random
	Random.seed!(n_wvfs_plot)
	n_wvfs = length(detector_plot_config_slider.efc_cut_left .< efc .< detector_plot_config_slider.efc_cut_right)
	wvf_index_efc_cut = collect(1:n_wvfs)[detector_plot_config_slider.efc_cut_left .< efc .< detector_plot_config_slider.efc_cut_right]
	n_wvfs = length(wvf_index_efc_cut)
	rand_wvfs_plot = rand(wvf_index_efc_cut, length(wvf_index_efc_cut))
end;

# ╔═╡ d8f70629-7019-48e6-b7cd-da5a0100bcda
begin
	@bind set_n Slider(1:n_wvfs - n_wvfs % n_wvfs_plot, show_value=true)
end

# ╔═╡ 77ddc928-4b62-4241-8001-01ff84969272
md" Selected subset: $set_n"

# ╔═╡ 54860a44-0823-4cb5-8958-474948138a25
begin
	idx_wvfs_plot = rand_wvfs_plot[set_n:set_n+n_wvfs_plot]
end;

# ╔═╡ 97f7da06-a687-4c62-a8d3-bc42430a9ff1
begin
	pe = plot(efc, st=:stephist, bins=detector_plot_config_slider.n_bins_efc, yscale=:log10, label="DAQ online energy", xlabel="Energy (ADC)", ylabel="Counts", size=(800, 500))
	plot!(title=format("{} Julia DSP Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
	ylims!(ylims()[1], ylims()[2])
	plot!(fill(detector_plot_config_slider.efc_cut_left, 2), [0.1, 1e4], label="Energy Cut Window", color=:red, lw=2.5, ls=:dot)
	plot!(fill(detector_plot_config_slider.efc_cut_right, 2), [0.1, 1e4], label="Energy Cut Window", color=:red, lw=2.5, ls=:dot, showlegend=false)
	p = plot(u"µs", NoUnits, size=(1000, 700), legend=:topright, thickness_scaling=1.0)
	xlims!(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right)
	# plot!(wvfs_fc[idx_wvfs_plot], label=permutedims(ts[idx_wvfs_plot]))
	for (iw, w) in enumerate(wvfs_type_selector)
		if iw == 1
			plot!(wvfs_options[w][idx_wvfs_plot], label=permutedims(ts[idx_wvfs_plot]), color=permutedims(wvfs_plot_colors), showlegend=true)
		else
			plot!(wvfs_options[w][idx_wvfs_plot], label=permutedims(ts[idx_wvfs_plot]), color=permutedims(wvfs_plot_colors), showlegend=false)
		end
	end
	for (iw, w) in enumerate(vline_type_selector)
		col = vline_plot_colors[iw]
		if occursin("area", w)
			for (area_max, wf, lab) in zip(vline_options[w][idx_wvfs_plot], wvfs_pz[idx_wvfs_plot], ts[idx_wvfs_plot])
				wf_area_trunc = TruncateFilter(area_max - dsp_config_slider.qdrift_window .. area_max)(wf)
				plot!(wf_area_trunc.time, fill(0.0, length(wf_area_trunc.signal)), fillrange=wf_area_trunc.signal, fillalpha=0.35, label=lab, color=col, showlegend=false)
			end
		end
		for (vline_type, lab) in zip(vline_options[w][idx_wvfs_plot], ts[idx_wvfs_plot])
			vline!([vline_type], label=lab, color=col, ls=:dot, lw=1.5, showlegend=false)
		end
	end
	plot!(title=format("{} Julia DSP Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
	PlutoPlot(plot([p, pe]..., layout=(1,2), size=(2500, 700), legend=:outertopright))
end

# ╔═╡ 12450361-8f9d-401b-9c44-0d9a4156251e
begin
	pt0 = stephist(t0, bins=0:0.1:120, label="t0", xlabel="t0", ylabel="Counts", xlims=(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right))
	pt50 = stephist(uconvert.(u"µs", t50), bins=0:0.1:120, label="t50", xlabel="t50", ylabel="Counts", xlims=(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right))
	pt80 = stephist(t80, bins=0:0.1:120, label="t80", xlabel="t80", ylabel="Counts", xlims=(detector_plot_config_slider.xlims_left, detector_plot_config_slider.xlims_right))

	plot([pt0, pt50, pt80]..., layout=(1,3), size=(2500, 700), legend=:outertopright, plot_title=format("{} Julia DSP Investigator ({}-{}-{}-{})", string(det), string(filekey.setup), string(filekey.period), string(filekey.run), string(filekey.category)))
end	

# ╔═╡ f52abfbf-1ddc-48c6-ab4d-406abf1069ee
begin
	ptrap = stephist(e_rtft, bins=0:detector_plot_config_slider.bin_width_e:maximum(e_rtft[isfinite.(e_rtft)]), xlabel="Energy (ADC)", ylabel="Counts", yscale=:log10, label="Trap", plot_title="Trap/ZAC/CUSP Optimized")
	stephist!(e_cusp, bins=0:detector_plot_config_slider.bin_width_e:maximum(e_cusp[isfinite.(e_cusp)]), xlabel="Energy (ADC)", ylabel="Counts", yscale=:log10, label="Cusp")
	stephist!(e_zac, bins=0:detector_plot_config_slider.bin_width_e:maximum(e_zac[isfinite.(e_zac)]), xlabel="Energy (ADC)", ylabel="Counts", yscale=:log10, label="ZAC")
	pe_qc = stephist(e_10410, bins=0:detector_plot_config_slider.bin_width_e:maximum(e_10410[isfinite.(e_10410)]), xlabel="Energy (ADC)", ylabel="Counts", yscale=:log10, label="Trap 10410", plot_title="10410/313/Max UnOptimized")
	stephist!(e_313, bins=0:detector_plot_config_slider.bin_width_e:maximum(e_313[isfinite.(e_313)]), xlabel="Energy (ADC)", ylabel="Counts", yscale=:log10, label="Trap 313")
	stephist!(wvf_max, bins=0:detector_plot_config_slider.bin_width_e:maximum(wvf_max[isfinite.(wvf_max)]), xlabel="Energy (ADC)", ylabel="Counts", yscale=:log10, label="Max")

	plot([ptrap, pe_qc]..., layout=(1,2), size=(2000, 700), legend=:outertopright)
end

# ╔═╡ Cell order:
# ╠═ecc1218c-a408-4de6-8cb5-40d29613b9e7
# ╟─98a808df-40ff-4eda-89ee-772147f7e42f
# ╟─f2cf0e3f-d544-4a0b-bcdf-d02bf9beb9d6
# ╟─2f337776-68cd-4b38-9047-88742bfa1c8a
# ╟─4933bfe7-2d34-4283-bd4f-2cb0006c5158
# ╟─d35a8320-ae1e-485d-8751-1a2b36f2b809
# ╟─24b387ff-5df5-45f9-8857-830bc6580857
# ╟─493bdebf-7802-4938-8693-5e0648bd8a2b
# ╟─a11bc8f3-c095-4100-96af-1e9d3eed0215
# ╟─8c632350-3fb5-462a-95f9-9b26306655fa
# ╟─39b468d4-673c-4834-b042-d80cd9e0dfe7
# ╟─2c4dd859-0aa6-4cd2-94b3-7a171368065b
# ╟─63b6be12-3b55-450a-8a2a-3d410f0dcb6c
# ╟─36232152-9811-4520-ba65-0f855e4ebbe4
# ╟─dd962859-e7a8-419e-87df-3c2499f013e5
# ╟─9edf7ca3-db5f-422d-ac11-f0e6fd17c087
# ╟─0e685d5c-2cfd-4f02-90a7-846b62a6426b
# ╟─5df81bc3-f86f-4a85-927c-ef9f3285dc97
# ╟─6f67faca-1065-43f6-94a2-345ac74a6a6f
# ╟─f5140e3a-2601-489a-9098-dc78b24ec0c3
# ╟─de7861dc-5665-4e88-ac00-af15ad1ada5c
# ╟─d8fef6df-d831-4248-aeae-84869714f76d
# ╟─73d945de-e5c7-427e-a75a-b0df28f86bd4
# ╟─4b7b64ec-9689-4a0a-acc1-8f90061e5498
# ╟─35fc04da-4adb-4af6-8f01-452b68577f87
# ╠═cf295ad2-13ca-4508-bf51-5cfc52229fbd
# ╟─73cbfa7e-b18c-4e2f-a39f-ae68bbf4a6fb
# ╟─e6a3a08d-9c64-451f-931c-54d8c3e3d6c5
# ╟─43204bc4-c869-4a76-ba93-500a33c39b89
# ╟─b5ddd0ba-a771-4f2e-94e1-ab9e994d69b3
# ╟─7ff8481d-4658-4a3f-a671-858edcf389cf
# ╟─f7211685-e046-47cc-8261-b0c54b466469
# ╟─56faaa8d-50e6-4b8e-aa7b-a8afd38e322d
# ╟─77ddc928-4b62-4241-8001-01ff84969272
# ╟─d8f70629-7019-48e6-b7cd-da5a0100bcda
# ╟─54860a44-0823-4cb5-8958-474948138a25
# ╟─97f7da06-a687-4c62-a8d3-bc42430a9ff1
# ╟─12450361-8f9d-401b-9c44-0d9a4156251e
# ╟─f52abfbf-1ddc-48c6-ab4d-406abf1069ee
# ╟─71f6a835-ebac-49e6-aac1-bbcc4d73bfe8
# ╟─89bfab7a-9f29-49e0-a313-544125367ee8
# ╟─ee95f35f-4c9c-4323-8f97-4beafab379fe
# ╟─bb81cf12-32c9-4bd0-92a8-b727fcf9b098
# ╟─60f190de-8ee1-4c81-a73d-6c5f1b52632e
# ╟─cc99755f-be83-4525-8850-a9df59eaca15
# ╟─d4d256f0-45c5-4f2b-bbc2-3f44b2802805
