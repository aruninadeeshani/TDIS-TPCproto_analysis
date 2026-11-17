#include "../TPC.C"
#include <chrono>
// root -l -b -q 'analysis_root/tests/ChannelCheck.C()'
// Focused unconnected/high-channel fixture using the installed ordered MS CSV.
void ChannelCheck() {
    using namespace PersonalTPC;
    auto folder=fs::temp_directory_path()/("tpc-channel-check-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder);
    struct Cleanup {fs::path path;~Cleanup(){fs::remove_all(path);}} cleanup{folder};
    auto write=[](const fs::path& path,const std::string& text){std::ofstream f(path);f<<text;};
    const int run=999101,largeRun=999102,observedRun=999103,sparseRun=999104;
    auto db=folder/"cuts.ini",study=folder/"study.ini",veto=folder/"veto.csv",binning=folder/"binning.csv",pairs=folder/"pairs.csv";
    write(db,"[global]\nsetup=medium_small\nreject=Time >= 228\n[run 999102]\nsetup=large_medium\n");
    write(study,"[global]\n");write(veto,"run,event\n999101,2\n");
    write(binning,"Branch,BinWidth,Bins,Min,Max,LogY\nChannel,1,,-0.5,,true\nPad,1,,-0.5,,true\n");
    write(pairs,"Name,X,Y,XBinWidth,XBins,XMin,XMax,YBinWidth,YBins,YMin,YMax,LogZ\n"
                "Time_vs_Channel,Channel,Time,,,,,,3,-0.5,299.5,false\n"
                "Pad_vs_FADC,FADC,Pad,,,,,,,,,false\n"
                "Custom_Channel,Channel,Time,,6,-0.5,5.5,,3,-0.5,299.5,false\n");
    Config config(folder.c_str(),(folder/"results").c_str(),veto.c_str(),db.c_str(),study.c_str());
    config.binningPath=binning;config.plots2DPath=pairs;
    const auto& ms=config.Mapping(run);
    SelfCheck(config,ms);
    Check(ms.pads.size()==506 && ms.blankTiles.size()==244,"Confirmed MS map is not installed");
    Check(ms.Number(576)==-1 && ms.Number(578)==-1 && ms.Number(289)==505 && ms.byChannel.rbegin()->first<=575,
          "Ordered MS CSV anchors or unconnected channels changed");
    Check(config.EffectiveChannels(run)==576 && config.EffectiveChannels(largeRun)==576,
          "Effective channel span did not follow the selected run map");
    Check(config.channelsPerModule==72 && config.readoutModules==8,"Configured module grouping changed");
    for(int channel:{576,578}) {
        Hit h{1,channel,1,10,100};
        Check(HitValue(h,"FADC",config,run)==channel+1 && HitValue(h,"Pad",config,run)==-1,
              "Unconnected channel lost its FADC alias or acquired a Pad ID");
        Check(PassHit(h,config,false,run) && PassHit(h,config,true,run),
              "Unconnected channel acquired an implicit hit cut");
    }
    auto aliases=folder/"aliases.ini";
    write(aliases,"[global]\nsetup=medium_small\nreject=FADC in 577,579\n");
    Config aliasConfig=config;aliasConfig.cutDBPath=aliases;aliasConfig.database.reset();aliasConfig.runInfo.clear();
    Check(!PassHit({1,576,1,10,100},aliasConfig,false,run) && !PassHit({1,578,1,10,100},aliasConfig,false,run) &&
          PassHit({1,575,1,10,100},aliasConfig,false,run),"FADC cuts are not one-based electronics aliases");
    write(aliases,"[global]\nsetup=medium_small\nreject=Pad in 0,505\n");
    aliasConfig.database.reset();aliasConfig.runInfo.clear();
    Check(!PassHit({1,152,1,10,100},aliasConfig,false,run) && !PassHit({1,289,1,10,100},aliasConfig,false,run) &&
          PassHit({1,576,1,10,100},aliasConfig,false,run) && PassHit({1,578,1,10,100},aliasConfig,false,run),
          "Mapped Pad cuts aliased unconnected channels");
    write(aliases,"[global]\nsetup=medium_small\nreject=Pad == -1\n");
    aliasConfig.database.reset();aliasConfig.runInfo.clear();
    Check(!PassHit({1,576,1,10,100},aliasConfig,false,run) && !PassHit({1,578,1,10,100},aliasConfig,false,run) &&
          PassHit({1,289,1,10,100},aliasConfig,false,run),"Explicit unmapped-Pad cut did not use Pad=-1");

    auto axis=RawBinning({152},"Channel",config,{},run);
    Check(axis.bins==576 && axis.low==-.5 && axis.high==575.5 && axis.logY,"Sparse MS input shrank the full Channel axis");
    axis=RawBinning({},"Channel",config,{},run);
    Check(axis.bins==576 && axis.high==575.5,"Empty MS input shrank the Channel axis");
    axis=RawBinning({152},"Channel",config,{},largeRun);
    Check(axis.bins==576 && axis.high==575.5,"LM axis inherited the MS span");
    axis=RawBinning({576,578},"Channel",config,{},run);
    Check(axis.bins==579 && axis.high==578.5,"Observed unconnected channels were clipped");
    axis=RawBinning({650},"Channel",config,{},run);
    Check(axis.bins==651 && axis.high==650.5,"Observed channels beyond the map were clipped");
    axis=RawBinning({578},"Channel",config,{{"Bins","6"},{"Min","-.5"},{"Max","5.5"}},run);
    Check(axis.bins==6 && axis.low==-.5 && axis.high==5.5,"Explicit Channel axis overrides were widened");
    auto custom=config;custom.binningPath=folder/"custom-binning.csv";
    write(custom.binningPath,"Branch,BinWidth,Bins,Min,Max,LogY\nChannel,1,6,-0.5,5.5,false\n");
    axis=RawBinning({578},"Channel",custom,{},run);
    Check(axis.bins==6 && axis.high==5.5 && !axis.logY,"Custom Channel CSV settings were overwritten");

    auto makeRun=[&](int number,const std::vector<Hit>& rows) {
        TFile file((folder/("ntuple"+std::to_string(number)+".root")).c_str(),"RECREATE");
        TTree tree("ntuple","Mapped and unconnected channel pulses");
        int event,channel,hit,time,charge;
        tree.Branch("Event",&event);tree.Branch("Channel",&channel);tree.Branch("Hit",&hit);tree.Branch("Time",&time);tree.Branch("Charge",&charge);
        for(const auto& h:rows){event=h.event;channel=h.channel;hit=h.hit;time=h.time;charge=h.charge;tree.Fill();}
        tree.Write();file.Close();
    };
    makeRun(run,{{1,576,1,10,100},{1,578,1,11,100},{1,576,2,228,100},{2,578,1,12,100},{1,289,1,13,100}});
    auto data=ReadRun(run,config);
    auto displayed=DisplayHits(data,1,config,ms,{});
    Check(displayed.size()==1 && displayed[0].channel==289,"Display included unconnected channels or dropped a mapped pad");
    for(bool cleaned:{false,true}) {
        Raw(config,run,cleaned,"Channel",&data);Raw(config,run,cleaned,"Pad",&data);Raw2D(config,run,cleaned,"",&data);
        auto out=config.OutputForRun(run)/"raw/run999101";if(cleaned)out/="cleaned";
        const int entries=cleaned?3:5,count=cleaned?1:2;
        TFile one((out/"branches_Channel/raw_histograms.root").c_str(),"READ");
        auto channel=dynamic_cast<TH1D*>(one.Get("raw_Channel"));
        Check(channel && channel->GetNbinsX()==579 && channel->GetEntries()==entries && channel->Integral()==entries,
              "Raw/cleaned observed Channel span or entry count failed");
        Check(channel->GetBinContent(channel->FindBin(576))==count && channel->GetBinContent(channel->FindBin(578))==count,
              "Raw/cleaned unconnected Channel bins lost hits");one.Close();
        TFile padFile((out/"branches_Pad/raw_histograms.root").c_str(),"READ");
        auto pad=dynamic_cast<TH1D*>(padFile.Get("raw_Pad"));
        Check(pad && pad->GetNbinsX()==506 && pad->GetEntries()==entries && pad->Integral()==1 &&
              pad->GetBinContent(0)==2*count && pad->GetBinContent(pad->FindBin(505))==1,
              "Raw/cleaned Pad bins lost unmapped underflow or aliased it to a connected pad");padFile.Close();
        TFile two((out/"2d/raw_2d_histograms.root").c_str(),"READ");
        auto channel2D=dynamic_cast<TH2D*>(two.Get("raw2d_Time_vs_Channel"));
        auto fadc=dynamic_cast<TH2D*>(two.Get("raw2d_Pad_vs_FADC"));
        auto narrow=dynamic_cast<TH2D*>(two.Get("raw2d_Custom_Channel"));
        Check(channel2D && channel2D->GetNbinsX()==579 && channel2D->Integral()==entries,
              "2D Channel axis dropped unconnected raw/cleaned hits");
        Check(fadc && fadc->GetBinContent(fadc->FindBin(577,ms.Number(576)))==count &&
              fadc->GetBinContent(fadc->FindBin(579,ms.Number(578)))==count && fadc->GetEntries()==entries &&
              fadc->Integral()==1 && fadc->GetBinContent(fadc->FindBin(290,505))==1,
              "2D FADC/Pad aliases lost unmapped underflow or a connected pad");
        Check(narrow && narrow->GetNbinsX()==6 && narrow->Integral()==0 &&
              narrow->Integral(7,7,0,narrow->GetNbinsY()+1)==entries,"2D custom Channel range did not retain overflow");two.Close();
    }
    auto summary=Analyze(run,100,config,ms);
    Check(summary.occupancy->GetNbinsX()==580 && summary.occupancy->GetXaxis()->GetXmax()==579.5,
          "Voltage occupancy ignored observed channels beyond the MS map");
    Check(summary.acceptedHits==4 && summary.time->Integral()==4 && summary.mappedHits==1 &&
          summary.pad->GetNbinsX()==506 && summary.pad->Integral()==1 && summary.vetoes.empty(),
          "Voltage analysis lost unconnected time hits or mapped them to pads");
    makeRun(sparseRun,{{1,152,1,10,100}});
    auto sparse=Analyze(sparseRun,100,config,ms);
    Check(sparse.occupancy->GetNbinsX()==577 && sparse.occupancy->GetXaxis()->GetXmax()==576.5,
          "Sparse voltage input shrank the selected MS channel span");
    std::vector<Hit> observed;
    for(int c=0;c<=650;++c)observed.push_back({1,c,1,10,100});
    makeRun(observedRun,observed);
    auto extended=Analyze(observedRun,100,config,ms);
    Check(extended.occupancy->GetNbinsX()==652 && extended.occupancy->Integral()==1 &&
          extended.occupancy->GetBinContent(extended.occupancy->FindBin(651))==1,
          "Voltage occupancy overflowed observed channels beyond the map");
    Check(extended.vetoes.size()==1 && extended.vetoes[0].detectorWide,"Extended-channel detector-wide threshold changed");

    // Exercise module 8 separately from module 0 and from adjacent events.
    RunData noise;noise.run=run;
    auto add=[&](int event,int channel,int time,int hit=1) {
        noise.events[event].push_back(noise.hits.size());noise.hits.push_back({event,channel,hit,time,100});
    };
    int highModule=8*config.channelsPerModule;
    for(int c=0;c<config.moduleChannels;++c) {
        add(1,highModule+c,40);add(2,c,40);
        add(3,(c%2?highModule:0)+c,40);
        add(4,highModule+c,40+c);
        for(int hit=1;hit<=config.moduleHitsPerChannel;++hit)add(5,highModule+c,1000*hit+c,hit);
    }
    add(6,1000000000,40);add(7,-1,40);
    auto burst=ClassifyEvent(noise,1,noise.events.at(1),config);
    Check(burst.modules.size()==1 && burst.modules[0].module==8 && burst.modules[0].channels==config.moduleChannels &&
          burst.modules[0].hits==config.moduleChannels && burst.modules[0].synchronous==config.moduleChannels,
          "Extended-module synchronous noise evidence changed");
    auto next=ClassifyEvent(noise,2,noise.events.at(2),config);
    Check(next.modules.size()==1 && next.modules[0].module==0 && next.modules[0].hits==config.moduleChannels,
          "Adjacent module-0 event collided with module 8");
    Check(!ClassifyEvent(noise,3,noise.events.at(3),config).Vetoed(),"Module 0 and module 8 occupancies were merged");
    Check(!ClassifyEvent(noise,4,noise.events.at(4),config).Vetoed(),"Candidate-only module became noise");
    auto repeated=ClassifyEvent(noise,5,noise.events.at(5),config);
    Check(repeated.modules.size()==1 && repeated.modules[0].module==8 && repeated.modules[0].synchronous==1 &&
          repeated.modules[0].hits==config.moduleChannels*config.moduleHitsPerChannel,"Extended-module repeated-hit threshold changed");
    Check(!ClassifyEvent(noise,6,noise.events.at(6),config).Vetoed(),"Sparse high module was rejected or aliased");
    bool failed=false;try{ClassifyEvent(noise,7,noise.events.at(7),config);}catch(const std::exception&){failed=true;}
    Check(failed,"Negative channel no longer rejected");

    // Exact adjacent-event regression: a flattened event*8+module key would
    // merge (event10,module8) with (event11,module0). Thresholds here are
    // local to the fixture; installed settings remain unchanged.
    Config three=config;three.moduleChannels=3;three.moduleSynchronous=3;
    RunData adjacent;adjacent.run=run;
    adjacent.hits={{10,576,1,40,100},{10,578,1,40,100},{11,0,1,40,100},{11,1,1,40,100}};
    adjacent.events[10]={0,1};adjacent.events[11]={2,3};
    Check(!ClassifyEvent(adjacent,10,adjacent.events.at(10),three).Vetoed() &&
          !ClassifyEvent(adjacent,11,adjacent.events.at(11),three).Vetoed(),
          "Adjacent events falsely merged sparse modules at threshold three");
    adjacent.events[10].push_back(adjacent.hits.size());adjacent.hits.push_back({10,577,1,40,100});
    auto positive=ClassifyEvent(adjacent,10,adjacent.events.at(10),three);
    Check(positive.modules.size()==1 && positive.modules[0].module==8 && positive.modules[0].channels==3 &&
          !ClassifyEvent(adjacent,11,adjacent.events.at(11),three).Vetoed(),
          "Three-channel module-8 burst affected adjacent module-0 event");
    std::cout<<"PASS: 506-pad MS map, unconnected Channels576/578, raw/cleaned 1D/2D and display exclusion, FADC/Pad aliases, custom axes, dynamic voltage occupancy and sparse module-noise grouping."<<std::endl;
}
