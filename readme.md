# **Patch**
Short spaced orthogonal lead investigation for patch-based ECG. The files in this repository contributed to several [Computing in Cardiology publications](https://scholar.google.com/citations?user=eHzaFo8AAAAJ&hl=en) as part of a PhD in Ulster University.

## **Structure**
### Data
The ECG data are from three sources:
1. Horacek data located in 'BalloonBSPMdata.mat'. This is 352 node baseline/annotation data
2. Kornreich BSPM data located in 'data_*.mat'. Three files with 352 node annotated beats
3. STAFF 9-lead ECGs (V1-V6;I-III) in '/STAFF/*.mat' with annotations in 'ann.xls'

Horacek and Kornreich are private, and can only be shared with permission from Ulster University. STAFF data is freely available on Physionet.

### Code
1. patch_v6.m - Master file
	- dataHorResponders.m - Sort patients by positive responders
	- steRank.m - Rank each possible bipolar lead by ST-elevation
	- filterLeadLength.m - Filter out leads over 100 mm
	- calculateSTE.m - Calculate the mean and median ST-elevation across remaining leads (all patients)
	- plotBSPM.m - Plot a BSPM contour map at the ST segment with SSLs annotated
	- steMaxIdx.m - Order patients by ST-elevation
	- plotSSL.m - Plot two leads on a subplot for a selected patient
	- splitData.m - Split data into training and test cell matrices
	- getCoeffs.m - Calculate linear regression coeffs for a single lead between two 352node BSPM points
	- getBspmSSL.m - Generate a single lead using linear regression coeffs and verify performance
	- splitDataAnn.m - Split data with annotations of MI, LVH, normal
	- compareLeads.m - Plot two leads to compare them
2. staffTest.m - Strip STAFF dataset to get annotations and individual beats
	- averageBeat.m - Returns 12-lead recording for one median beat
	- detectSTEMI.m - Returns STEMI classification from median 12-lead beat and their J-points
	- qrsDetector.m - Returns R-wave locations 
3. horacekFormatting.m - Returns 12-lead ECGs from the horacek dataset


### Use
Please note, the STAFF III database must be downloaded before use. The ann.xlsx file must be in the same directory.
