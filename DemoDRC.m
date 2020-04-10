function [M,M2,P,MortR] = DemoDRC(Amin)
% Returns contact matrix for community and home absed on specified
% compression
% Input
% Amin - The minimum age of the different classes
% Output
%M - contact matrix (Size: AxA)
%M2 - contact matrix home (Size: AxA)
%P -Population size


MortRFull=[0.0026 
    0.0026
0.0148
0.0148
0.0600
0.0600
0.146
0.146
0.295
0.295
1.25
1.25
3.99
3.99
8.61
8.61
13.4
13.4
13.4
13.4
13.4]./100;

% Double checked pfull
PFull=[15502762
13298331
11123785
9110315
7460648
6229691
5163631
4271737
3506681
2878188
2339248
1856309
1430823
1064142
763636
460778
224804
81333
20241
3316
168];
AA=5.*[0:(length(PFull)-1)];
Findx=cell(length(Amin),1);
F80indx=cell(length(Amin),1);
for ii=1:length(Amin)-1
    gg=find(AA==Amin(ii));
    hh=find(AA==Amin(ii+1));
    Findx{ii}=[gg:(hh-1)];
    F80indx{ii}=[gg:(hh-1)];
end
Findx{ii+1}=hh:length(PFull);
gg=find(AA==Amin(end));
hh=find(AA==80);
F80indx{ii+1}=[gg:(hh-1)];

% Double checked MFULL (All locations)
MFull=[5.216865469	2.717833541	1.38868468	0.779944047	1.132602294	1.63305069	1.725589778	1.335339743	0.694045134	0.389796193	0.379386276	0.307503165	0.191652424	0.137708668	0.065036696	0.031544293
2.499477349	14.00086465	2.863396677	0.786885533	0.50899388	1.130658394	1.354400371	1.268936722	0.884352474	0.400682763	0.272567124	0.23839157	0.181984743	0.113744104	0.046896869	0.027299546
0.940158667	4.310282469	12.88011987	1.575155419	0.757351165	0.697119715	0.823768961	0.964941313	0.874492943	0.479498216	0.266356359	0.142059893	0.097409042	0.093279002	0.056124992	0.02981955
0.537236553	1.170793696	4.355679437	10.19076275	2.102505022	1.035892566	0.697435501	0.857380903	0.765445102	0.637463458	0.350355411	0.149642013	0.093252434	0.056570546	0.027618914	0.014485231
0.836134677	0.736464356	0.800720008	3.678777937	5.047920172	2.370081724	1.406364734	1.044553262	0.784479319	0.787079984	0.538079672	0.257727923	0.115091454	0.038787611	0.029824557	0.019118914
1.352712303	0.831613455	0.45000759	1.138036891	2.72265641	3.535718443	2.054393973	1.438061275	1.047935934	0.820535189	0.712344055	0.30276605	0.156579509	0.044368557	0.01530591	0.01235305
1.428929045	2.058085005	1.424749785	0.65079132	1.320972952	2.174208949	2.515635596	1.860748751	1.293198871	0.960004147	0.677723431	0.361731178	0.203217794	0.055916513	0.025982866	0.016321934
1.196883845	1.816483651	1.329787639	0.749740675	0.775687839	1.443198322	1.743822614	2.092641999	1.594527255	1.011470255	0.717061429	0.272408138	0.173881607	0.081814526	0.036553606	0.011841468
0.838392372	1.358107721	1.31746288	1.01588554	0.817172614	1.103892507	1.430908802	1.479360852	1.593950703	1.171999063	0.824661132	0.219615761	0.181632526	0.072461586	0.03192433	0.014584378
0.527131091	1.023161114	0.949156869	1.135191925	0.706055457	0.874224288	1.075837352	1.168040094	1.085635012	1.075961294	0.75075513	0.280529616	0.128262588	0.043158836	0.027763602	0.023497456
0.527997918	0.914159248	1.03137581	0.867529754	0.697776232	0.980426527	1.000368931	0.938794895	1.146468513	1.170731758	0.907498192	0.394829864	0.150815817	0.037428134	0.027199632	0.021624339
1.05284688	1.393980507	0.876688672	0.734120318	0.536430164	0.830206986	0.858051723	0.619212026	0.608288081	0.589532302	0.578503349	0.316672367	0.206549085	0.091113968	0.028957496	0.03065894
0.664085652	0.757724389	0.49376257	0.420049781	0.317844496	0.445982157	0.508028462	0.55526258	0.437881126	0.383347687	0.348620186	0.248083835	0.144799218	0.092917065	0.034762984	0.009647471
0.510279601	0.83647296	0.67495253	0.258969819	0.235856263	0.257521319	0.354015782	0.366174127	0.273039616	0.152768147	0.170710301	0.166931327	0.126743826	0.107168164	0.063750809	0.017041927
0.212997884	0.681720867	0.582453864	0.387669485	0.102610005	0.194619527	0.150458727	0.261294556	0.235086558	0.217337445	0.20520921	0.130967156	0.139382812	0.101691476	0.072826419	0.063572241
0.291516805	0.429453097	0.57319755	0.402306512	0.097546024	0.108397824	0.129515774	0.212211262	0.187870035	0.188999316	0.204912678	0.13275535	0.058920351	0.087386023	0.049473263	0.032118124
];
% Double checked MFULL2 (HOME)
M2Full=[0.529288913	0.976502695	0.661260874	0.292391633	0.41056708	0.554563178	0.53988024	0.485754522	0.26098452	0.138510616	0.145311057	0.139114386	0.090538211	0.045974451	0.013892348	0.012527289
0.571141585	0.773627219	0.793862126	0.323820297	0.141937078	0.373722647	0.475783478	0.460045242	0.343406822	0.178577066	0.122735706	0.132749338	0.095979873	0.047402644	0.019006366	0.011221572
0.458418709	0.939139339	1.148996307	0.507290603	0.159919572	0.163809033	0.246876461	0.360231576	0.33611086	0.201008052	0.113298297	0.064167516	0.052721317	0.047744479	0.027589394	0.01022935
0.294004248	0.499903169	0.755502253	0.779087408	0.311055724	0.167726561	0.105725497	0.256510579	0.26566002	0.260635286	0.179239277	0.079523663	0.056551032	0.032626734	0.015321989	0.008032437
0.52011932	0.381723074	0.345760391	0.597044292	0.580158317	0.435802864	0.159693084	0.089966622	0.146464826	0.232144016	0.173808585	0.111681452	0.043274576	0.017663111	0.010103456	0.00695123
0.841655808	0.437220654	0.220154616	0.251783445	0.403016705	0.454447113	0.261982213	0.10065701	0.051248392	0.099953964	0.155475233	0.099034773	0.056381904	0.023128552	0.004331265	0.006879486
0.936542417	1.168185043	0.739641	0.180507014	0.245416027	0.422750227	0.314736273	0.263134397	0.134471414	0.079054949	0.091135478	0.1036443	0.083821046	0.024917238	0.012183538	0.006183377
0.793523003	1.12678315	0.904876685	0.379814202	0.110796542	0.1366873	0.200371308	0.287925061	0.15228907	0.071457373	0.064061032	0.050234115	0.067444317	0.03961127	0.015394819	0.005217977
0.640539708	1.0034482	0.934809501	0.518671356	0.188548809	0.122549534	0.175245075	0.209306584	0.173341348	0.125214845	0.068158031	0.029442216	0.074484923	0.051187455	0.017723137	0.009833942
0.371190283	0.697766395	0.734150212	0.575416416	0.2993499	0.157968998	0.095793397	0.155591971	0.123337205	0.147183827	0.12687423	0.060739928	0.045232882	0.028716434	0.017558667	0.01726734
0.445860332	0.447686635	0.562007258	0.379265994	0.281669774	0.208969591	0.132031942	0.08853079	0.089594937	0.113226031	0.110657912	0.107183864	0.047228479	0.022948923	0.018623396	0.017048821
0.904038602	1.03124619	0.594805273	0.457777233	0.317210465	0.382321583	0.338206654	0.16366329	0.103037095	0.204903006	0.202083279	0.119514728	0.12020364	0.073539942	0.019427573	0.026136787
0.580790595	0.666858789	0.427990288	0.276743517	0.146828192	0.162011724	0.182461352	0.172405735	0.111149179	0.075092088	0.108476485	0.104513143	0.061582033	0.058191708	0.019822129	0.004677847
0.46600921	0.757470872	0.637125114	0.23316182	0.165957644	0.155389787	0.224600953	0.273425177	0.203880055	0.104840666	0.133360294	0.118929615	0.088937211	0.061843477	0.039369179	0.010786598
0.198721294	0.655742442	0.54924355	0.324136124	0.055793338	0.130262383	0.089666662	0.182935369	0.163693548	0.159042123	0.145217832	0.088110535	0.098259944	0.06137425	0.03816691	0.034458307
0.274099618	0.39313431	0.553266064	0.371152489	0.077873312	0.076056111	0.086458135	0.183264293	0.166718143	0.168672642	0.192811929	0.122671792	0.049873633	0.074422663	0.039340152	0.026682198];

P=zeros(length(Amin),1);
Ptemp=zeros(length(MFull(:,1)),length(Amin));
PtempMort=zeros(length(MortRFull(:,1)),length(Amin));

Mtemp=zeros(length(Amin),length(MFull(:,1)));
M2temp=zeros(length(Amin),length(MFull(:,1)));

for ii=1:length(P)
    P(ii)=sum(PFull([Findx{ii}]));    
    Ptemp([F80indx{ii}],ii)=PFull([F80indx{ii}])./sum(PFull([F80indx{ii}]));
    PtempMort([Findx{ii}],ii)=PFull([Findx{ii}])./sum(PFull([Findx{ii}]));
    Mtemp(ii,:)=sum(MFull([F80indx{ii}],:),1);
    M2temp(ii,:)=sum(M2Full([F80indx{ii}],:),1);
end
MortR=(MortRFull')*PtempMort;
M=Mtemp*Ptemp;
M2=M2temp*Ptemp;

M=(M+M')./2;
M2=(M2+M2')./2;

end
