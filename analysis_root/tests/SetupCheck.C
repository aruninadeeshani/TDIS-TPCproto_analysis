#include "../TPC.C"
#include <chrono>
// root -l -b -q 'analysis_root/tests/SetupCheck.C()'
// All MS electronics below are SYNTHETIC, not a proposed wiring map.
void SetupCheck() {
    using namespace PersonalTPC;
    auto folder=fs::temp_directory_path()/("tpc-setup-check-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder/"source/auxiliary");
    struct Cleanup {fs::path path;~Cleanup(){fs::remove_all(path);}} cleanup{folder};
    auto write=[](const fs::path& path,const std::string& text){std::ofstream f(path);f<<text;};
    auto contents=[](const fs::path& path){std::ifstream f(path);return std::string(std::istreambuf_iterator<char>(f),{});};
    auto dbPath=folder/"cuts.ini",study=folder/"study.ini",veto=folder/"veto.csv",binning=folder/"binning.csv",pairs=folder/"pairs.csv";
    write(study,"[global]\n");write(veto,"run,event\n");
    write(binning,"Branch,BinWidth,Bins,Min,Max,LogY\nPad,1,,-0.5,,true\nTime,1,,0,100,false\nEvent,1,1,0,50000,false\n");
    write(pairs,"Name,X,Y,XBinWidth,XBins,XMin,XMax,YBinWidth,YBins,YMin,YMax,LogZ\nTime_vs_Pad,Pad,Time,,,,,,10,0,100,false\n");
    write(dbPath,"[global]\n");
    Config defaults("","","",dbPath.c_str());
    Check(defaults.Info(1615).setup=="large_medium" && defaults.Info(1616).setup=="large_medium","Historical runs changed default setup");
    Check(!defaults.Info(1616).driftVoltageV && !defaults.Info(1616).tiltAngleDeg,"Unknown metadata invented");
    Check(defaults.Mapping(1616).pads.size()==410 && defaults.Mapping(1616).Number(504)==0,"Default LM mapping changed");
    fs::copy_file(defaults.Mapping().path,folder/"source/auxiliary/pad_map.csv");

    write(dbPath,"[run 01616]\nsetup = medium_small\ndrift_voltage_v =\ntilt_angle_deg = 0\n"
                 "[global]\nsetup = large_medium\ndrift_voltage_v = -1.2e3\ntilt_angle_deg = -2.5\n"
                 "[run 1617]\nsetup = medium_small\ntilt_angle_deg =\n");
    CutDatabase inherited(dbPath);
    Check(inherited.Info(1615).setup=="large_medium" && *inherited.Info(1615).driftVoltageV==-1200,"Global inheritance failed");
    Check(inherited.Info(1616).setup=="medium_small" && !inherited.Info(1616).driftVoltageV && *inherited.Info(1616).tiltAngleDeg==0,"Run override/empty reset failed");
    Check(*inherited.Info(1617).driftVoltageV==-1200 && !inherited.Info(1617).tiltAngleDeg,"Independent metadata inheritance failed");
    for(const auto& bad:std::vector<std::string>{
        "[global]\nsetup=medium_small\nsetup=large_medium\n",
        "[global]\ndrift_voltage_v=\nDrift_Voltage_V=1\n",
        "[run 1616]\ntilt_angle_deg=0\ntilt_angle_deg=\n",
        "[global]\nsetpu=large_medium\n", "[global]\nsetup=\n", "[global]\nsetup=small\n",
        "[global]\ndrift_voltage_v=NaN\n", "[global]\ntilt_angle_deg=inf\n",
        "[global]\ndrift_voltage_v=1e999\n", "[global]\ntilt_angle_deg=12 deg\n",
        "[run 01616]\n[run 1616]\n"}) {
        write(dbPath,bad);bool failed=false;
        try{CutDatabase invalid(dbPath);}catch(const std::exception& e){failed=std::string(e.what()).find(dbPath.string()+":")!=std::string::npos;}
        Check(failed,"Invalid metadata accepted or missing line diagnostic");
    }

    // Real LM medium channels are reused at a different Pad origin; small
    // channels are arbitrary unused electronics, assigned only for this test.
    std::set<int> used;std::ostringstream map;
    map<<"Pad,FADC,Channel,PadSize,Row,Col,X_mm,Y_mm,Width_mm,Height_mm\n";
    for(int i=0;i<250;++i) {
        int channel=defaults.Mapping().pads.at(160+i).channel;used.insert(channel);
        map<<i<<','<<channel+1<<','<<channel<<",Medium,"<<i/25<<','<<i%25<<','<<4*(i%25)<<','<<4*(i/25)<<",4,4\n";
    }
    int nextChannel=0,pad=250;
    for(int row=0;row<10;++row)for(int col=12;col<39;++col) {
        if(row==0 && col<24)continue; // Twelve interior/edge blanks.
        while(used.count(nextChannel))++nextChannel;
        used.insert(nextChannel);
        map<<pad++<<','<<nextChannel+1<<','<<nextChannel<<",Small,"<<row<<','<<col<<','<<2*col<<','<<40+2*row<<",2,2\n";
        ++nextChannel;
    }
    Check(pad==508,"Bad synthetic fixture count");
    auto synthetic=folder/"source/auxiliary/pad_map_medium_small.csv";
    write(synthetic,map.str());
    write(dbPath,"[global]\nsetup=large_medium\ndrift_voltage_v=1200\ntilt_angle_deg=-3\nreject=Pad == 0\n"
                 "[run 1616]\nsetup=medium_small\ndrift_voltage_v=\ntilt_angle_deg=0\n");
    Config config(folder.c_str(),"",veto.c_str(),dbPath.c_str(),study.c_str());
    config.source=folder/"source";config.output=folder/"results";config.binningPath=binning;config.plots2DPath=pairs;
    const auto& lm=config.Mapping(1615);const auto& ms=config.Mapping(1616);
    Check(lm.Number(152)==160 && ms.Number(152)==0,"Same Channel did not change Pad with setup");
    Check(!PassHit({1,504,1,10,10},config,false,1615) && PassHit({1,152,1,10,10},config,false,1615),"LM Pad0 alias changed");
    Check(!PassHit({1,152,1,10,10},config,false,1616) && PassHit({1,504,1,10,10},config,false,1616),"MS Pad0 cut used LM alias");
    Check(ms.pads.size()==508 && ms.width==100 && ms.height==60 && ms.blankTiles.size()==242,"MS connected/blank geometry wrong");
    Check(ms.pads[250].section=="Small" && ms.pads[507].number==507,"Connected small numbering wrong");
    std::vector<int> coverage(6000);
    auto tiles=ms.pads;tiles.insert(tiles.end(),ms.blankTiles.begin(),ms.blankTiles.end());
    for(const auto& p:tiles) {
        if(p.number<0)Check(p.channel==-1 && p.fadc==-1 && p.section=="Small","Blank acquired an electronics/Pad ID");
        for(int y=int(p.y);y<p.y+p.height;++y)for(int x=int(p.x);x<p.x+p.width;++x)++coverage.at(y*100+x);
    }
    Check(std::all_of(coverage.begin(),coverage.end(),[](int n){return n==1;}),"Small physical footprint has gaps/overlaps");
    auto axis=RawBinning({0},"Pad",config,{},1616);
    Check(axis.bins==508 && axis.low==-.5 && axis.high==507.5 && axis.logY,"Automatic Pad extent depends on observed hits");
    axis=RawBinning({-1,507},"Pad",config,{{"Bins","5"},{"Min","-.5"},{"Max","9.5"}},1616);
    Check(axis.bins==5 && axis.high==9.5,"Explicit axis overrides ignored");
    Check(config.OutputForRun(1615)==config.output && config.OutputForRun(1616)==config.output/"setups/medium_small","Default setup output contract changed");
    Check(config.OutputForRuns({1615,1616})==config.output/"setups/mixed","Mixed output not isolated");
    Config overrideLM(folder.c_str(),"",veto.c_str(),dbPath.c_str(),study.c_str(),"large_medium");overrideLM.output=config.output;
    Check(overrideLM.Info(1616).setup=="large_medium" && overrideLM.Mapping(1616).Number(504)==0 && overrideLM.OutputForRun(1616)==config.output,"Manual LM override failed");
    Config explicitMS(folder.c_str(),(folder/"explicit").c_str(),veto.c_str(),dbPath.c_str(),study.c_str(),"medium_small");
    Check(explicitMS.OutputForRun(1615)==folder/"explicit" && explicitMS.OutputForRuns({1615,1616})==folder/"explicit","Explicit output changed");
    Config missing(folder.c_str(),"",veto.c_str(),dbPath.c_str(),study.c_str(),"medium_small");
    missing.source=folder/"absent";bool failed=false;
    try{missing.Mapping(1615);}catch(const std::exception& e){failed=std::string(e.what()).find("pad_map_medium_small.csv")!=std::string::npos;}
    Check(failed,"Missing confirmed MS map silently fell back");

    int unconnected=0;while(ms.Number(unconnected)>=0)++unconnected;
    for(int run:{1615,1616}) {
        TFile file((folder/("ntuple"+std::to_string(run)+".root")).c_str(),"RECREATE");
        TTree tree("ntuple","Synthetic setup pulses");
        int event=run==1615?42535:5584,channel,hit=1,time=10,charge=100;
        tree.Branch("Event",&event);tree.Branch("Channel",&channel);tree.Branch("Hit",&hit);tree.Branch("Time",&time);tree.Branch("Charge",&charge);
        for(int c:{152,ms.pads.back().channel,unconnected}){channel=c;tree.Fill();}
        tree.Write();file.Close();
    }
    for(bool cleaned:{false,true}) {
        Raw(config,1616,cleaned,"Pad");Raw2D(config,1616,cleaned);
        auto out=config.OutputForRun(1616)/"raw/run1616";if(cleaned)out/="cleaned";
        TFile one((out/"branches_Pad/raw_histograms.root").c_str(),"READ");
        auto h=dynamic_cast<TH1D*>(one.Get("raw_Pad"));
        Check(h && h->GetNbinsX()==508 && h->GetEntries()==(cleaned?2:3),"Raw/clean Pad count or selection wrong");
        Check(h->GetBinContent(1)==(cleaned?0:1) && h->GetBinContent(508)==1 && h->GetBinContent(0)==1,"Raw/clean mapping bins wrong");
        auto setup=dynamic_cast<TNamed*>(one.Get("setup")),drift=dynamic_cast<TNamed*>(one.Get("drift_voltage_v"));
        Check(setup && std::string(setup->GetTitle())=="medium_small" && drift && std::string(drift->GetTitle()).empty(),"UNCUT metadata provenance missing");
        Check(contents(out/"branches_Pad/pad_map_used.csv")==map.str(),"1D snapshot is not effective map");
        one.Close();
        TFile two((out/"2d/raw_2d_histograms.root").c_str(),"READ");
        auto h2=dynamic_cast<TH2D*>(two.Get("raw2d_Time_vs_Pad"));
        Check(h2 && h2->GetNbinsX()==508 && h2->GetEntries()==(cleaned?2:3),"2D setup/bin/cleaning mismatch");two.Close();
        Check(contents(out/"2d/pad_map_used.csv")==map.str(),"2D snapshot is not effective map");
    }
    auto lmSummary=Analyze(1615,100,config,ms),msSummary=Analyze(1616,100,config,lm);
    Check(lmSummary.pad->GetNbinsX()==410 && msSummary.pad->GetNbinsX()==508 && msSummary.timePad->GetNbinsX()==508,"Voltage histograms share a fixed map");
    lmSummary.selected=msSummary.selected=true;failed=false;
    try{Overlay({lmSummary,msSummary},"bad_overlay",false,true,config);}catch(...){failed=true;}
    Check(failed,"Cross-setup Pad overlay accepted");
    Examples(config);
    auto exampleDir=config.output/"setups/mixed/comparisons";
    Check(contents(exampleDir/"run1616_pad_map_used.csv")==map.str(),"Examples reused first run map");
    Check(fs::exists(config.output/"events/run1615/event_1615_42535.root") &&
          fs::exists(config.output/"setups/medium_small/events/run1616/event_1616_5584.root"),"Examples output paths are not per run");
    write(dbPath,"[global]\nsetup=large_medium\ndrift_voltage_v=-100.25\nreject=Pad == 0\n"
                 "[run 1616]\nsetup=medium_small\n");
    Config scan=config;scan.database.reset();scan.runInfo.clear();
    Check(ReadVoltages(scan).at(1616)==-100.25,"Voltage did not inherit signed fractional DB value");
    Voltage(scan);
    for(int run:{1615,1616}) {
        auto out=config.output/"setups/mixed/voltage"/config.Info(run).setup;
        TFile f((out/"time_histograms.root").c_str(),"READ");
        auto h=dynamic_cast<TH1D*>(f.Get(("pad_run"+std::to_string(run)).c_str()));
        Check(h && h->GetNbinsX()==(run==1615?410:508),"Mixed voltage winner/map lost");
        auto info=dynamic_cast<TNamed*>(f.Get(("run"+std::to_string(run)+"_drift_voltage_v").c_str()));
        Check(info && std::string(info->GetTitle())=="-100.25","Voltage export truncated/lost signed DB value");f.Close();
    }
    // Public entry point accepts the final setup parameter; requesting MS
    // without a source map is tested above, while LM stays publicly runnable.
    TPC("info",1616,-1,folder.c_str(),(folder/"public").c_str(),veto.c_str(),dbPath.c_str(),study.c_str(),"large_medium");
    std::cout<<"PASS: native metadata validation/inheritance, overrides, Pad0 aliases, dynamic bins, blank geometry, raw/clean provenance, mixed examples and voltage separation."<<std::endl;
}
