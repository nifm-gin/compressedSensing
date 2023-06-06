% Script for making undersampling patterns --------------------------------

%% Add BART v0.6.00 path
% set Matlab path and TOOLBOX_PATH environment variable
%bartPath = '/opt/bart-0.6.00';
bartPath = '~/Apps/Bart/bart-0.8.00';
addpath(fullfile(bartPath, 'matlab'));
setenv('TOOLBOX_PATH', bartPath);
bart('version')

addpath(pwd,'diffusionKPCA/KernelLib' );
addpath(pwd, 'diffusionKPCA/ReconLib' );
addpath(pwd, 'diffusionKPCA/PreImLib' );
addpath(pwd, 'diffusionKPCA/FFTLib' );
addpath(pwd, 'diffusionKPCA' );

%% 2023.03.15 - Undersampling patterns for Tonton RV 
close all
clear
clc

dataInfo.AF = [2,4];  %  acceleration factor (possibilité de generer avec plusieurs valeurs séparés par une virgule)
dataInfo.marginAF = -0.0
; % facteur a retirer au nombre de ligne total pour que ça ne plante pas...)
dataInfo.undersamplingPatternType = ["MonteCarlo-VariableDensity","Poisson-UniformDensity-EllipticalScanning"]; % type de mask : "MonteCarlo-VariableDensity","Poisson-UniformDensity-EllipticalScanning","Poisson-VariableDensity-EllipticalScanning","Poisson-UniformDensity","Poisson-VariableDensity"];dataInfo.undersamplingPatternType = ["MonteCarlo-VariableDensity","Poisson-UniformDensity-EllipticalScanning","Poisson-VariableDensity-EllipticalScanning","Poisson-UniformDensity","Poisson-VariableDensity"];%"MonteCarlo-VariableDensity","Poisson-UniformDensity-EllipticalScanning","Poisson-VariableDensity-EllipticalScanning","Poisson-UniformDensity","Poisson-VariableDensity"]; (possibilité de generer avec plusieurs valeurs séparés par une virgule)
dataInfo.undersamplingPatternStrategy = ["Multi3D"]; % Multi3D, Mono3D (possibilité de generer avec plusieurs valeurs séparés par une virgule)
%dataInfo.undersamplingPatternStrategy = ["Multi3D"]; % Multi3D, Mono3D (possibilité de generer avec plusieurs valeurs séparés par une virgule)
dataInfo.undersamplingPatternDimensions = [360,180,220,4,7]; % Readout x Phase1 x Slice selection x Number of channels x Number of volumes (direction Bo + Diff_Dir)
%dataInfo.undersamplingPatternDimensions = [180,90,110,4,1]; % Anat = 1 seul volume.
%dataInfo.undersamplingPatternDimensions = [360,180,220,4,1]; % Anat = 1 seul volume.
dataInfo.undersamplingPatternFullySampledRegionArea = [(pi*(15)^2)]; % surface du disque FullySampled.
dataInfo.B0SamplingStrategy = ["UndersampledB0"]; % FullySampledB0 or UndersampledB0 = obligatorie pour l'Anat (possibilité de generer avec plusieurs valeurs séparés par une virgule) 
dataInfo.B0Number = 1; % 1, car 1 voulume b=0, (pour l'anat,  = 0)

dataInfo.showResults = 1; % 1 pour visualiser les resultat du calcul du mask avant de sauvegarder, sinon 0 
dataInfo.slicesToShow = [90,45,55,1,1]; % on a ffiche les coupes centrales seulement


saveUndersamplingPatterns =1; % 0 just calculate and display - 1 also write Mask files
pathToSaveUndersamplingPatterns = "/home_ldap/_SHARE/Compressed-Sensing/Codes/CS_DTI/Masks/";





%%

makeSaveUndersamplingPatterns_v5(dataInfo,pathToSaveUndersamplingPatterns,saveUndersamplingPatterns);

