clear all;
close all;

% get(0,'DefaultFigureWindowStyle')
set(0,'DefaultFigureWindowStyle','docked')
% set(0,'DefaultFigureWindowStyle' , 'normal')

addpath(genpath('D:\Data\Chronux\'))

probeNr = [0];
Area= {'V1';'LM'};
sampleDef = [];
reps=10;

for i = 1:length(probeNr)
    pathName = 'X:\data\Ai9xCD1_68\Ai9xCD1_68_20230124_Recording1\9_FGinverseBlue_g0\';
    paradigmName = '1_stimulusStageGratingFigureInverse_rgb';

    load([pathName paradigmName '\mergedCombined_Uncurated_DataForProbe' num2str(probeNr(i)) '.mat']);
    rS = getConditionParameterMatrixForParadigm('stimulusStageGratingFigureInverse_rgb', params);
    %     sampleDef.R = rS.R;
    showAllParadigms(params);

    stimInfo = load([pathName paradigmName '\stimulation.mat']);

    edges = [(-499):1:2500];

    try
        stims.sampleInds = stims.sampleIndsStart;
    end

    stims.sampleInds = round(stims.sampleInds/30);

    lfpdata = LFP_1kHz;

%     if i==1
%         [SpikesV1, qualitiesV1, depthsV1] = loadSpikeDataQuality(SpikeTimes*30000,SpikeIds,SpikeQualities, SpikeDepths, uniqueAbsoluteSpikeInds,30);
%         wfV1=wf;
%         uniqueAbsoluteSpikeIndsV1=uniqueAbsoluteSpikeInds;
%     elseif i==2
%         [SpikesLM, qualitiesLM, depthsLM] = loadSpikeDataQuality(SpikeTimes*30000,SpikeIds,SpikeQualities, SpikeDepths, uniqueAbsoluteSpikeInds,30);
%         wfLM=wf;
%         uniqueAbsoluteSpikeIndsLM=uniqueAbsoluteSpikeInds;
%     end

    dirb=rS.dirBackground;
    dirf=rS.dirFigure;
    radii=rS.R;
    FigContrast=rS.contrastFigure;
    BgContrast=rS.contrastBackground;
    RedAmps=rS.RedlaserAmpl;
    BlueAmps=rS.BluelaserAmpl;
    radii1=[20 45]/2;
    radii2=[10 20 30 45 60]/2;
    radii3=[10 20 30 45 60]/2;

    for db=1:length(dirb)
        for df=1:length(dirf)
            for r1=1:length(radii1)
                for la=1:length(BlueAmps)

                sampleDef.dirBackground=dirb(db);
                sampleDef.dirFigure=dirf(df);
                sampleDef.contrastBackground = 110;
                sampleDef.contrastFigure = 110;
                sampleDef.Stim0Inhib1 = 1;
                sampleDef.RedlaserAmpl = 0;
                sampleDef.BluelaserAmpl = BlueAmps(la);
                sampleDef.R=radii1(r1);

%                 if i==1
%                     [spikedattrialsV1,sV1,V1inds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesV1,edges);
%                     V1FigGrdSpikes(:,:,:,db,df,r1) = spikedattrialsV1(:,:,1:reps);%time,units,reps,db,df,la
%                     V1indsmat(:,db,df,r1) = V1inds(1:reps);
%                 elseif i==2
%                     [spikedattrialsLM,sLM,LMinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesLM,edges);
%                     LMFigGrdSpikes(:,:,:,db,df,r1) = spikedattrialsLM(:,:,1:reps);%time,units,reps,db,df,la
%                     LMindsmat(:,db,df,r1) = LMinds(1:reps);
%                 end

                [FigGrdLFPdattrials,FigGrndsLFP,FigGrndLFPTrialinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,lfpdata,edges);
                FigGrdLFP(i,:,:,:,db,df,r1,la) = FigGrdLFPdattrials(:,:,1:reps);%probe,time,chn,reps
                FigGrdTrialindsmat(:,db,df,r1,la) = FigGrndLFPTrialinds(1:reps);
                end
            end
        end
    end

    clear  sampleDef
    dirf2=[0:pi/2:(2*pi-pi/2)];
    for df2=1:length(dirf2)
        for r2=1:length(radii2)
            for la=1:length(BlueAmps)

            sampleDef.dirBackground=pi;
            sampleDef.dirFigure=dirf2(df2);
            sampleDef.contrastBackground = 0;
            sampleDef.contrastFigure = 110;
            sampleDef.Stim0Inhib1 = 1;
            sampleDef.RedlaserAmpl = 0;
            sampleDef.BluelaserAmpl = BlueAmps(la);
            sampleDef.R=radii2(r2);

%             if i==1
%                 [FigspikedattrialsV1,sV1,V1Figinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesV1,edges);
%                 V1FigSpikes(:,:,:,df2,r2) = FigspikedattrialsV1(:,:,1:reps);%time,units,reps,db,df,la
%                 V1Figindsmat(:,df2,r2) = V1Figinds(1:reps);
%             elseif i==2
%                 [FigspikedattrialsLM,sLM,LMFiginds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesLM,edges);
%                 LMFigSpikes(:,:,:,df2,r2) = FigspikedattrialsLM(:,:,1:reps);%time,units,reps,db,df,la
%                 LMFigindsmat(:,df2,r2) = LMFiginds(1:reps);
%             end

            [FigLFPdattrials,FigsLFP,FigLFPTrialinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,lfpdata,edges);
            FigLFP(i,:,:,:,df2,r2,la) = FigLFPdattrials(:,:,1:reps);%probe,time,chn,reps
            FigTrialindsmat(:,df2,r2,la) = FigLFPTrialinds(1:reps);
            end
        end
    end

    clear  sampleDef
    dirb2=[pi,pi/2];
    for db2=1:length(dirb2)
        for r3=1:length(radii3)
            for la=1:length(BlueAmps)

            sampleDef.dirBackground=dirb(db2);
            sampleDef.dirFigure=0;
            sampleDef.contrastBackground = 110;
            sampleDef.contrastFigure = 0;
            sampleDef.Stim0Inhib1 = 1;
            sampleDef.RedlaserAmpl = 0;
            sampleDef.BluelaserAmpl = BlueAmps(la);
            sampleDef.R=radii3(r3);

%             if i==1
%                 [HolespikedattrialsV1,sV1,V1Holeinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesV1,edges);
%                 V1HoleSpikes(:,:,:,db2,r3) = HolespikedattrialsV1(:,:,1:reps);%time,units,reps,db,df,la
%                 V1Holeindsmat(:,db,df) = V1Holeinds(1:reps);
%             elseif i==2
%                 [HolespikedattrialsLM,sLM,LMHoleinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesLM,edges);
%                 LMHoleSpikes(:,:,:,db2,r3) = HolespikedattrialsLM(:,:,1:reps);%time,units,reps,db,df,la
%                 LMHoleindsmat(:,db2,r3) = LMHoleinds(1:reps);
%             end

            [HoleLFPdattrials,HolesLFP,HoleLFPTrialinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,lfpdata,edges);
            HoleLFP(i,:,:,:,db2,r3,la) = HoleLFPdattrials(:,:,1:reps);%probe,time,chn,reps
            HoleTrialindsmat(:,db2,r3,la) = HoleLFPTrialinds(1:reps);
            end
        end
    end
    clear  sampleDef
    for la=1:length(BlueAmps)

    sampleDef.contrastBackground = 0;
    sampleDef.contrastFigure = 0;
    sampleDef.Stim0Inhib1 = 1;
    sampleDef.RedlaserAmpl = 0;
    sampleDef.BluelaserAmpl = BlueAmps(la);

%     if i==1
%         [V1Con0dattrials,V1sCon0,V1Con0Trialinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesV1,edges);
%         V1Con0Spikes(:,:,:) = V1Con0dattrials(:,:,1:reps);%time,units,reps,laser
%         V1Con0Trialindsmat(:) = V1Con0Trialinds(1:reps);
%     elseif i==2
%         [LMCon0dattrials,LMsCon0,LMCon0Trialinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,SpikesLM,edges);
%         LMCon0Spikes(:,:,:) = LMCon0dattrials(:,:,1:reps);%time,units,reps,laser
%         LMCon0Trialindsmat(:) = LMCon0Trialinds(1:reps);
%     end

    [Con0LFPdattrials,Con0sLFP,Con0LFPTrialinds] = getParadigmTrials('stimulusStageGratingFigureInverse_rgb',stimInfo,stims,sampleDef,lfpdata,edges);
    Con0LFP(i,:,:,:,la) = Con0LFPdattrials(:,:,1:reps);%probe,time,chn,reps
    Con0Trialindsmat(:,la) = Con0LFPTrialinds(1:reps);
    end
    clear  sampleDef
end