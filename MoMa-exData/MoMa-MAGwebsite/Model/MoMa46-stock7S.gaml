/**
 *  MoMa 4.6 
 *  Author: MANEERAT Somsakun, UMR IDEES, University of Rouen
 *
 *  Description: MoMa (Model Of Mosquito Aedes aegypti), an agent-based model of Aedes aegypti female mosquitoes,
*   provides spatially explicit information on mosquito behaviour and aims to lead new research questions and target actions against dengue mosquitoes.
 *  This model is a thesis project developped in the framwork of French national program ANR AEDESS (http://anr-aedess.fr) 
 *  and European program FP7 DENFREE (http://www.denfree.eu/)
 *
  *  UPDATE: 05/05/2015 
 */

model MoMa46_STOCK
global
{

	/**********************
	 * CHECKING PARAM
	 ***********************/
	 bool bOpenMole <- false parameter: "use the model with OpenMole?" category 'Sim Options';
	 bool bAvgDenseBSMode <- true parameter: "use average density for BS?" category: 'Sim Options'; //depend on landuse_class.csv (surveyHzrnSaket.csv) = do not use min &max for SS
	 bool bAutoSimStop <- false parameter: "auto stopping?" category: 'Sim Options';
	 bool bCheckMaxDisp <- true parameter: "check max dispersal?"  category: 'Sim Options';	
	 bool bCheckStck <- true parameter: "check Aedes'stocks?"  category: 'Sim Options';
	 bool bUpdateStck <- true parameter: "update aquatic stocks?"  category: 'Sim Options';
	 bool bCheckActStat <- false parameter: "check activites?" category: 'Sim Options';
	 bool bCheckNbBiteOvipoPerHourInObj <- false parameter: "Check nb bites/ovipos per hous in each spatial obj?" category: 'Sim Options';
	 bool bCheckFinalFuzzyScore <- false parameter: "show final fuzzy score?" category: 'Sim Options';
	 bool bDisprsShpNeed <- false parameter: "need dispersal shapefile/ show dispersal circle?"  category: 'Sim Options';
 	 bool bNeedNewZone <- true parameter: "need new study zone shapefile?" category: 'Sim Options';
	 bool bShowLimitArea <- false parameter: "show limit area for Aedes?" category: 'Sim Options';
	 bool bTrackAedes <- false parameter: "track Aedes movement?" category: 'Sim Options';
	 bool bTrackFuzzyProcess <- false parameter: "Track fuzzy process? " category: 'Sim Options';
	 bool bVarComputeExecTime <- false parameter: "Compute and display the execution time of each reflex?" category: 'Sim Options';
	
	/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	 * OPENMOLE Parameters for SpatObj (Scenario 1 - test P1/P4/P5) only
	 *xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
	float fgRecCapaRate  <- 0.1 parameter: "max human density? " category: 'Variable to test'; //P1 - blood (value: 0.1 - 9.1, step 0.5)
	float fgDenBsInNoWater <- 0.0 parameter: "density BS inside? " category: 'Variable to test'; //P4 - BS In (value: 0.0 - 2, step 0.05)
	float fgDenBsOutNoWater <- 0.0 parameter: "density BS outside? " category: 'Variable to test'; //PS - BS Out (value: 0.0 - 2, step 0.5)
		

	/********************************
	 * SIM ENVIRONMENT DATA
	 ********************************/
	 //spatial data
	 file SpatObj_shape <- file ("../shp/HzrSaket/ClassifHzrnSaket-survey-pop.shp") parameter: "Shapefile for the study zone" category: "GIS" ;
	 file SpatObjType_file <- csv_file("../shp/HzrSaket/surveyHzrnSaket.csv") parameter: "land use classes for the shapefile" category: "GIS";
	 geometry shape <- envelope(SpatObj_shape);
	 
	 
	//-------------- meteo data-------------------------	
	file temperatureFile <- csv_file("../num_data/meteo/meteo2008.csv") parameter: "temperature file" category: "Weather";
	file sunFile <- csv_file("../num_data/meteo/sun2008.csv") parameter: "sunset-sunrise file" category: "Weather";
	matrix temperature_matrix <- matrix(temperatureFile);
	matrix sun_matrix <- matrix(sunFile) ;
	matrix temperature_nbdays_dev_matrix <- matrix(csv_file("../num_data/meteo/EFFECT_TEMP_DEVDURATION.csv",";")) ;
	matrix<float> avgHourAirTempK_matrix; // matrix to memorise the average air temperature for each hour of the total days in simulation

	//------------- name of output files -> to modify while changing the landuse csv file------------
	string dispersShpName <- "dispersal46-OM_S1.shp"; 
	string MaxDisPopShpName <- "Maxdispersal46-OM_S1.shp";
	string newMapName <- "simHzrnsaket46-OM_S1.shp";
	string deathFileName <- "death46_OMS1.csv";
	string fileDailyStock <- "OsDailyStcks46_OMS1.csv";
	string fileEggLaid <- "OsEggLaid46_OMS1.csv";
	string fileFuzzyProc <- "fuzzyTarget46-OM_S1.csv";
	string fileSelBestTarget <- "bestTarget46-OM_S1.csv";

	string txtSaveDeath <- "";
	string txtStock <- "";
	string txtEggLaid <- "";
	
	/****************
	 * SCENARIOS 
	 ****************/
	 int NB_ADULT_FEMALE <- 100 parameter: "Nb female per a breeding site" category: "Aedes scenario"; 
	 list<int> listOrigBs <- [0,3954,2428,1203,3352,1345,2079] parameter: "ID of Breeding site lists" category: "Aedes scenario"; 
	 
	/*************************************************
	 * 				AEDES'PARAMETERS
	 *************************************************/
	 float seed <- 354.0 parameter: true;
 	
	 //----- actions'time control ---------
	 float MAX_NEEDTIME_NECTAR <- 60.0; //s
	 float MIN_NEEDTIME_NECTAR <- 30.0; 
	 float MAX_NEEDTIME_BLOOD <- 60.0; //1mn
	 float MIN_NEEDTIME_BLOOD <- 5.0; //5sec
	 float MIN_NEEDTIME_LAYEGG <- 33.0*60; //mn*sec
	 float MAX_NEEDTIME_LAYEGG <- 120.0*60; //mn*sec	 
	 
	 //------- energy control -------------
	float NECTAR_GAIN_ENERGY_SPEED <- 1080.0/60; //1080 /60 (mn of energy/sec)
	float BLOOD_GAIN_ENERGY_SPEED <- (86400*3)/60; //= 4320sec = 86400sec/day * 3 days of survival for 1 full meal
	
	float TAKENECTAR_ENERGLOST_RATE <- 0.9; // Almeida, 2011
	float TAKEBLOOD_ENERGLOST_RATE <- 0.9; // almeida, 2011
	float WANDER_ENERGLOST_RATE <- 0.6; //almeida, 2011

	float ENERGY_LOST_LAYEGGE_SEC <- 1080.0; // to verified
	float AVG_CONSUM_ENERGY_DAY <- 0.60; // average of energy consumption per day = 60% (Templin, 2000; Almeida, 2011)
	float MAX_ENERGY <- 86400* 3* AVG_CONSUM_ENERGY_DAY; // seconds = the number of days that a mosquito can survive without feeding 3 days
	float BLOOD_GONOUSE_ULPERHOUR <- 0.3; //0.3ul/hour (Christopher, 1960 - full meal (3ul) is eliminate within 10hours in general) 
	float BLOOD_USE_LAYEGG_RATE <- 0.03; // estimate 100 eggs / 3mg (full blood meal)
	

	 //------ takeBlood probability ---------
	 float BITE_ATTEMP_SURV <- 0.98; // 2% Almeida p.1503
	 float BITING_SUCC_RATE <- 0.99; 
	 float BITING_SPEED_UL_SEC <- 0.05; //ul/sec : 1ul = 1mg
	 float LOWERLIMIT_MIN_BLOOD_NEED_DEVEGG <- 0.5; //mg the variable will use to calculate fMinBloodNeedDevEgg in Aedes
	 float UPPERLIMIT_MIN_BLOOD_NEED_DEVEGG <- 1.0; //mg the variable will use to calculate fMinBloodNeedDevEgg in Aedes
	 float AVG_MIN_BLOOD_NEED_DEVEGG <- 0.8; //mg the variable will use to calculate fMinBloodNeedDevEgg in Aedes
	 int MIN_NBBITES_GONO <- 1;
	 int MAX_NBBITES_GONO <- 9;	
	 
	 //------ lay egg probability -----------
	 float MAX_SPEED_LAYEGG <- 0.05; //0.05eggs/s = 3eggs/mn
	 float MIN_SPEED_LAYEGG <- 0.016; //0.016 eggs/s = 1egg/mn = 33eggs mn/oviposit
	 float MIN_TEMP_LAYEGG <- 18.0;//Â°C
	 //------- re-bite probability -----------
	 float MAX_STOCK_BLOOD <- 3.0;//maximum 3 mg/1 gono cycle
	 
	//-------- flight probability --------------
	float TARGET_MAXDIST_TOLERENCE <- 1.0;//0.1; //10cm
	float MAX_FLIGHT_SPEED <- 1.0; //meter/sec
	float MIN_FLIGHT_SPEED <- 0.5; //meter/sec
	float PERCEPTION_RADIUS <- 10.0; //meters
	//------- age -----------
	float MAX_AGE_ADULT <- 30.0; //30
	float MIN_AGE_ADULT <- 21.0; //21 A VOIR 
	float MATING_PROB <- 0.95;
	float DIST_NOEFF <- 0.1; //= 10cm distance which not effect the target decision
	
	//-------- physiologic dev -----
	float RO25_GONO <- 0.00898 ;
	float DHA_GONO <- 15725.23 ;
	float DHH_GONO <- 1756481.0 ;
	float THALF_GONO <- 447.17 ;	
	float R <- 1.987 ;
	float FRS_GONO_THRES <- 1.0 ;
	float LATE_GONO_THRES <- 0.58 ;

	float hourlyPhysioDevGono;
		
	//-------- Stages'development ------------
	float NOMINAL_DAILYSURV_EGGS_RATE <- 0.99;
	float NOMINAL_DAILYSURV_LARVAE_RATE <- 	0.99;
	float NOMINAL_DAILYSURV_PUPAE_RATE	<- 0.99;
	float NOMINAL_DAILYSURV_FEMALE_ADULTS_RATE	<- 0.91; // 0.091% for mortality in Otero, 2010
	
	// the value comes from Enzyme kenetics equation use to initiat an aquatic stock list in SpatObj
	int MAX_DAYNEED_EMBRYONATION <- 8; //days
	int MAX_DAYNEED_PUPATION <-	30; //days
	int MAX_DAYNEED_EMERGENCE	<-7	; //days
	int MIN_DAYNEED_EMBRYONATION<-	2; //days
	int MIN_DAYNEED_PUPATION	<-3	; //days
	int MIN_DAYNEED_EMERGENCE	<-1	; //days
	
	//------------ Egg hatching --------------
	float EGG_FLOOD_HATCH_RATIO_MAX <- 0.72; // value from CimSim
	float EGG_FLOOD_HATCH_RATIO_MIN <- 0.48; //value from CimSim
	float EGGS_HATCH_NO_FLOODING_PROB  <- 0.197; //Skeeter buster
	
	//--------- stage transition survival rate ----
	//--- values from Almeida - 2010 (who cited the value from Bellows et al. 1992,  Hawkins et al, 1997, Mwangi and Rembold, 1988, Peter and Barbosa, 1977)
	float EMBRYONATION_SURV	<- 0.7; //0.7
	float PUPATION_SURV	<- 0.7; //0.7;
	float EMERGC_SURV <- 0.7; //0.7;
	float FAMALE_PROB <- 0.5; //0.5;//prob to become a female adult
	
	//--------- survival rate with KTemp
	float EGG_DAILYSURV_LOW_TEMP_RATE <- 0.05 ;
	float EGG_DAILYSURV_HIGH_TEMP_RATE <- 0.05 ;
	
	float LARVAE_DAILYSURV_LOW_TEMP_RATE  <- 0.05 ;
	float LARVAE_DAILYSURV_HIGH_TEMP_RATE <- 0.05 ;

	float PUPAE_DAILYSURV_LOW_TEMP_RATE  <- 0.05 ;
	float PUPAE_DAILYSURV_HIGH_TEMP_RATE <- 0.05 ;
	
	float ADULT_DAILYSURV_LOW_TEMP_RATE <- 0.05 ;
	float ADULT_DAILYSURV_HIGH_TEMP_RATE  <- 0.05 ;
	
	//---------- CTemp limit for different state survival
	float EGG_DAILYSURV_LOW_TEMP_LIMIT <- -14.0 ;
	float EGG_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT  <- -6.0 ;
	float EGG_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT <- 30.0 ;
	float EGG_DAILY_SURVIVAL_HIGH_TEMP_LIMIT <- 47.0 ;
	
	float LARVAE_DAILYSURV_LOW_TEMP_LIMIT <- 5.0;
	float LARVAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT  <- 10.0;
	float LARVAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT  <- 39.0;
	float LARVAE_DAILYSURV_HIGH_TEMP_LIMIT <- 44.0;	
	
	
	float PUPAE_DAILYSURV_LOW_TEMP_LIMIT <- 5.0;
	float PUPAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT <- 10.0;
	float PUPAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT  <- 39.0;	
	float PUPAE_DAILYSURV_HIGH_TEMP_LIMIT <- 44.0;	

	float ADULT_DAILYSURV_LOW_TEMP_LIMIT <- 0.0 ;
	float ADULT_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT <- 4.0 ;
	float ADULT_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT  <- 40.0 ;
	float ADULT_DAILYSURV_HIGH_TEMP_LIMIT <- 50.0 ;		
	
	
	//paramerter for choice equation
	float BLOODEFF_WG_4_BLOOD <- 0.7;
	float BLOODEFF_WG_4_NECTAR <- 0.1; //dep more to energy level
	float BLOODEFF_WG_4_SHADE <- 0.5;
	float BLOODEFF_WG_4_BS <- 0.9;
	
	/////////// LAND USE //////////////
	float MAX_WORLDTEMPC <- 50.0;
	float AVG_EGG_BSSURFACE <- 300.0; //max 750eggs/container avg size (Iquito) Wong J, Stoddard ST, Astete H, Morrison AC, et al. (2011) Oviposition Site Selection by the Dengue Vector Aedes aegypti and Its Implications for Dengue Control. PLoS Negl Trop Dis 5(4): e1015. doi:10.1371/journal.pntd.0001015; http://www.plosntd.org/article/info:doi/10.1371/journal.pntd.0001015
	//AVG_EGG par m² a verifier la valeur car c'est a peu pres 70eggs par m2 (la valeur sortir de modèle) entre le mois may et oct dans Hopp, M.J and Foley J.A ( Global-scale relationships between climate and the dengue fever vector, aedes aegypti
	
	float MIN_TEMP_USECOOLER <- 20.0;//A CHANGER
	
	
	/**********************************************************
	 *             SIMULATION VARIABLES
	 **********************************************************/

	//------------- time variables --------------
	int nb_minutes update: int(time/ 60); // /60s = gama's variable! It represents the current simulated time in seconds
	int nb_hours update : int (nb_minutes/60); // 
	int nb_days update: int(nb_hours/ 24); 
	int hour_of_a_day update: nb_hours mod 24; // count an hour for each day (in pacific standard time  0 to 24h)
	int mn_of_an_hour update: nb_minutes mod 60; // count a minute for each hour
	bool bNewDay function: {hour_of_a_day mod 24 = 0 and mn_of_an_hour mod 60 = 0 and cycle != 0? true: false};
	bool bNewHour function: {mn_of_an_hour mod 60 = 0 and time mod 60 = 0 and cycle != 0? true: false};
	float curTime function: {hour_of_a_day + (mn_of_an_hour/100)};
	int sim_nb_daysInYear <- 154; //1/06/2008
	int iSimRealDay update: nb_days + sim_nb_daysInYear;
	int sim_days_duration <- 60; // How many days lates the simulation?
	int iNextDaySunriseHour;
	string strCurrentDay;
	
	
	int totalDataDays;	
	list<float> list_daily_max_temp <- [];
	list<float> list_daily_min_temp <- [];
	list<float> list_daily_max_wind <- [];
	list<float> list_daily_rain <- [];
	list<float> list_sunset  <- [];
	list<float> list_sunrise <- [];
	float prevAirTempHour <- 999.9; // temperature per hour		
	float fAirHourAvgTempK ;
	float fGlobalTodayAirTempC;
	float fGlobalTodayRainfall;

	map<string, int> mapStateIndex; //memorize the keyword and the position in other lists
	map<string, int>  mapTargetIndex;	
	list<string> list_KeyTargetName;
	
	list<float> list_percentPresPopResid; //percentage of human presences on residential surface by hour 
	list<float> list_percentPresPopComm; //percent of human presences on public (commercial) surface by hour 
	list<float> list_percentPresPopTrans; // percentage of human presences in a transport surface by hour
	list<float> list_AeActifRateByHour <- [0.05,0.05,0.05,0.2,0.6,0.95,0.95,0.85,0.7,0.6,0.3,0.1,0.05,0.1,0.3,0.5,0.65,0.85,0.95,0.95,0.6,0.2,0.1,0.05]; //percentage of aedes actif by hour, calculate based on the sunrise at 5 am
	
	float time_rest;
	float time_stat_takeBlood;
	float time_nectar;
	float time_lay_egg;
	
	float step <- 1°mn; // how many seconds per step for the simulation?; // 1 step = 1mn (60s) 
			
	int mn_sunrise;
	int hr_sunrise;
	
	list<float> stat_takeBlood <- [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
	list<float> stat_nectar <- [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
	list<float> stat_rest <- [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
	list<float> stat_wander <- [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
	list<float> stat_flyTo <-[0.0,0.0,0.0,0.0,0.0,0.0,0.0];
	
	int iNbAvgDailyTb;
	int iNbAvgDailyTn;
	int iNbAvgDailyBs;
	int iNbAvgDailyShd;
	int iNbAvgDailyFt;
	int iNbAvgDailyWd;
	float fAvgDailyFlyDist;	
	
	//save activity for Aedes
	map activities_color <- ["nectar" :: rgb("green"), "blood" :: rgb("red"), "shade" :: rgb("magenta"), "bs"::rgb("blue"), "flyTo"::rgb("yellow"),"wander"::rgb("orange")];
	bool debug <- false;
	 
	float totEggs;
	float totLarvae;
	float totPupae;
	int totNewAdults;
					
	float fMaxDist;
	
	float t1;
	float t2;
	float t3;
	float t4;
	float t5;
	float t6;
	float t7;
	float t8;
	float t9;
	float t10;
	float t11;
	float t12;
	float t13;
	float t14;
	float t15;
	float t16;
	
	/************************
	 * initialized processes
	 ************************/
	 
	 
	init
	{
		mapStateIndex <- ["virgin"::0, "ovipo"::1, "gono"::2];
		list_KeyTargetName <- ["nectar", "blood", "shade", "bs"];
		mapTargetIndex <- [list_KeyTargetName[0]::0, list_KeyTargetName[1]::1, list_KeyTargetName[2]::2,list_KeyTargetName[3]::3];
		
		list_percentPresPopResid <- [0.94,0.94,0.94,0.9,0.8,0.7,0.7,0.6,0.5,0.3,0.2,0.2,0.3,0.2,0.2,0.2,0.3,0.35,0.4,0.55,0.7,0.7,0.8,0.9];
		list_percentPresPopComm <- [0.01,0.01,0.01,0.01,0.05,0.1,0.1,0.1,0.2,0.5,0.7,0.7,0.55,0.7,0.75,0.75,0.6,0.5,0.4,0.3,0.2,0.2,0.1,0.05];//% par rapport au nb pop init pour les zones public
		list_percentPresPopTrans <-[0.05,0.05,0.05,0.09,0.15,0.20,0.20,0.30,0.30,0.20,0.10,0.10,0.15,0.10,0.05,0.05,0.10,0.15,0.20,0.15,0.10,0.10,0.10,0.05];
		
		do loadDailyWeather; //initialize totalDataDays at the same time
		avgHourAirTempK_matrix <- 0.0 as_matrix({totalDataDays,24});
		do createHourAirTempMatrix;
		
		hourlyPhysioDevGono <- 0.0;		
		fMaxDist <- 0.0;
		
		mn_sunrise <- int((list_sunrise at (nb_days) - int(list_sunrise at (nb_days)))*100);
		hr_sunrise <- int(list_sunrise at (nb_days));		
		
		time <- (hr_sunrise*60 + mn_sunrise)* 60.0; //initiate the sim time to the sunrise of the first day in seconds	
		iNextDaySunriseHour <- int(list_sunrise at (nb_days+1));
		fGlobalTodayAirTempC <- (list_daily_max_temp[nb_days] + list_daily_min_temp[nb_days])/2;
		fGlobalTodayRainfall <- list_daily_rain[nb_days];	
		 
		do createSpatObj;
		
		//create scenarios where Aedes born in different BS
		if(listOrigBs != [])
		{
			loop i from : 0 to: length(listOrigBs)-1
			{
				do createAedes (listOrigBs[i]);				
			}
		}
			
		if(bCheckActStat)
		{
			create Stat;
		}
			
	}
	
	/*--------------------------------------------------------
	*checking computational time in each function
	-------------------------------------------------------*/
	reflex info when: bVarComputeExecTime {
		if (t15 = 0) {t15 <- machine_time;}
		write "**********" + cycle + "*********";
		write "total time :" +(machine_time - t15);
		write "temps compute_duration_stat: " + t1;
		write "update_stock: " + t2;
		write "updateCone : " + t3;
		write "updateGonoPhysioDev : " + t4;
		write "checkCurrentAct : " + t5;
		write "update_limit : " + t6;
		write "saveMosAct : " + t7;
		write "updateDynmTargetsAttrc : " + t9;
		write "updateBloodAttract : " + t10;
		write "checkAedesActifPeriod : " + t11;
		write "action fuzzy	: " + t16;
		
		write "updateStocksAedes: " + t12;
		
		write "create_act_geom: " + t13;
		write "updateControlParam : " + t14;
	}
	
	
	/*****************************************
	 * STOP Simulation if bAutoSimStop = true
	 *****************************************/
	reflex stopSimulation when: bAutoSimStop and length(Aedes) = 0 //(nb_days = sim_days_duration or int(total_duration) > 72000000 //nb_days = sim_days_duration or length(Aedes) = 0
	{
		if(bCheckActStat and length(activity_stat) != 0)
		{
			save activity_stat type:"shp" to: "actTask.shp" with:[aedesName::"AEDES", type::"TYPE", duration::"DURATION", beginning::"START", sim_day::"DAY_SIM", hour::"HOUR", mn::"MN"];
		}
		
		if(bDisprsShpNeed and length(Dispersal) != 0)
		{
			save Dispersal type:"shp" to: MaxDisPopShpName with:[shape::shape, fMaxRadius::"MAXDIST"];
		}
		
		if(bTrackAedes and length(dailyTrail) != 0)
		{
			save dailyTrail type:"shp" to: "dailyTrail.shp" with:[shape::shape, aedesName::"AEDES"];			
		}
		save txtSaveDeath type: "csv" to: deathFileName;	
		save txtStock type: "csv" to: fileDailyStock;
		save txtEggLaid type: "csv" to: fileEggLaid;
		//save [average_duration, total_duration] type: "csv" to: "sim_time.csv";
		do halt;
	}	
	
	reflex mosqNominalProbKiller when: bNewDay
	{
		//write "nb days " + nb_days;
		ask Aedes
		{
			 if(!flip(NOMINAL_DAILYSURV_FEMALE_ADULTS_RATE))
			 {
			 	do killMe("random");//iAge will be update after in agent Aedes
			 }
		}
	}

		
	action createSpatObj
	{

		do loadClassLanduse;
		create SpatObj from: SpatObj_shape with: [strClass::string(read("LAND_USE"))]
		{
			OBJID <- name;
			spcClass <- (Class_Landuse ) first_with (each.strClassName = strClass);
			if (length(shape.geometries) > 1) 
			{
				loop i from: 1 to:  length(shape.geometries) - 1
				{
					create SpatObj {
						strClass <- myself.strClass;
						spcClass <- myself.spcClass;
						shape <-myself.shape.geometries[i];
						shape <- (shape + 0.05) simplification 0.01;
						OBJID <- name;
					}
				}
				shape <- shape.geometries[0];
			}
			shape <- (shape + 0.05) simplification 0.01;
		}
		
		
		// ----- correction the holes in each polygon------
		ask SpatObj 
		{
			if (shape != nil and shape.area > 0) 
			{
				do initialisation;
			} else {
				do die;
			}
		}
		//----- save new shape file for afterward analysis -----
		if(bNeedNewZone)
		{
			save SpatObj type:"shp" to: newMapName with:[OBJID::"OBJID", strClass::"strClass"];
		}
	
	}
	
	/*****************************************************************
	 *  import characteristics of different land uses
	 *********************************************************/
	action loadClassLanduse
	{
		
		
		if(bAvgDenseBSMode) //si oui, create land use class avec la valeur moyenne du bs in and out
		{
			create Class_Landuse from: SpatObjType_file header: true with:
			[
				strClassName::string(read("strClassName")),
				fRecCapaRate::float(read("fRecCapaRate")),
				socialClass::string(read("socialClass")),
				fsurfPrivRate::float(read("fsurfPrivRate")),
				fSurfPubRate::float(read("fsurfPubRate")),
				fSurfTransRate::float(read("fSurfTransRate")),
				fNectarSourceRate::float(read("fNectarSourceRate")),
				fRestSourceRate::float(read("fRestSourceRate")),
				fOutSpaceRate::float(read("fOutSpaceRate")),
				fDenBsOutNoWater::float(read("fDenBsOutNoWater")),
				fDenBsInNoWater::float(read("fDenBsInNoWater")),
				color::rgb(read("color")),
				fTempVarRate::float(read("fTempVarRate")),
				fPorosityMin::float(read("fPorosityMin")),
				fPorosityMax::float(read("fPorosityMin"))
			];
		}
		else
		{
			create Class_Landuse from: SpatObjType_file header: true with:
			[
				strClassName::string(read("strClassName")),
				fRecCapaRate::float(read("fRecCapaRate")),
				socialClass::string(read("socialClass")),
				fsurfPrivRate::float(read("fsurfPrivRate")),
				fSurfPubRate::float(read("fsurfPubRate")),
				fSurfTransRate::float(read("fSurfTransRate")),
				fNectarSourceRate::float(read("fNectarSourceRate")),
				fRestSourceRate::float(read("fRestSourceRate")),
				fOutSpaceRate::float(read("fOutSpaceRate")),
				fDenBsOutNoWaterMin::float(read("fDenBsOutNoWaterMin")),
				fDenBsOutNoWaterMax::float(read("fDenBsOutNoWaterMax")),
				fDenBsInNoWaterMin::float(read("fDenBsInNoWaterMin")),
				fDenBsInNoWaterMax::float(read("fDenBsInNoWaterMax")),
				color::rgb(read("color")),
				fTempVarRate::float(read("fTempVarRate")),
				fPorosityMin::float(read("fPorosityMin")),
				fPorosityMax::float(read("fPorosityMin"))
			];
		}
		
		if(bOpenMole) // si c'est pour le test avec OpenMole, prend la valeur de variable globale
		{
			ask Class_Landuse
			{
				fRecCapaRate <- fgRecCapaRate;
				fDenBsOutNoWater <- fgDenBsOutNoWater;
				fDenBsInNoWater <- fgDenBsInNoWater;
				//fPorosityAvg <- fgPorosityAvg;
			}
		}
	}
	
	
	/**********************************************
	 * adjust display color for different land use
	 **********************************************/
	rgb adjustObjColor(rgb myCol)
  	{
  		rgb newCol <- rgb("red");
  		switch (myCol)
  		{
  			match rgb("magenta") {newCol <- rgb(255,204,204);}
  			match rgb("yellow") {newCol <- rgb(255,229,204);}
  			match rgb("red") {newCol <- rgb(91,12,12);}
  			match rgb("brown") {newCol <- rgb(204,102,0);}
  			match rgb("green") {newCol <- rgb(229,255,204);}
  			match rgb("grey") {newCol <- rgb(224,224,224);}
  		}
  		return newCol;
  	}
  		
  	
  	/*************************************************
  	 * instantiate aedes agents in different scenarios
  	 *************************************************/	
	action createAedes(int numBs)
	{
		create Aedes number: NB_ADULT_FEMALE
		{
			location <- SpatObj[numBs].location;
			bornObj <- (SpatObj ) first_with (each.location = self.location);
			bornLoc <- location;
			perception_area <- shape + PERCEPTION_RADIUS;
			iAge <- 0; 
		}
		
		//create Dispersal shape if required
		if(bDisprsShpNeed)
		{
			ask SpatObj[numBs] 
			{
				create Dispersal 
				{
					strOrigSpatObj <- SpatObj[numBs].name;
					myself.myMaxDisp <- self;
				}
			}
		}
	}	
		
	/*****************************************
	 * load daily weather data from CSV file
	 * ******************************************/
	action loadDailyWeather
	{
		totalDataDays <- length(temperature_matrix column_at 0) - 1; //minus the title line
		loop i from: sim_nb_daysInYear to: (totalDataDays)
		{
			add (float((temperature_matrix column_at 1) at i)) to: list_daily_max_temp;	//c
			add (float((temperature_matrix column_at 2) at i)) to: list_daily_min_temp;	//c
			add (float((temperature_matrix column_at 3) at i)) to: list_daily_max_wind;	//km/h
			add (float((temperature_matrix column_at 4) at i)) to: list_daily_rain;	// mm
			add (float((sun_matrix column_at 1) at i)) to: list_sunrise;	
			add (float((sun_matrix column_at 2) at i)) to: list_sunset;	
		}
	}
	

	/****************************************************************
	 * calculate the average air temperature for each hour of a day
	 ***************************************************************/
	action createHourAirTempMatrix
	{
		loop i from: sim_nb_daysInYear to: (totalDataDays - 1)//a line  = a day
		{
			loop j from: 0 to: 23 // column for each hour
			{
				avgHourAirTempK_matrix[i,j] <- calHourlyTempK(j,"air");
			}
		}	
	}
	

	/***************************************************************************
	 * Calculate hourly air & water temperature by using TM model (for air temperatur only)
	 **************************************************************************/
	float calHourlyTempK (float calTime, string tempFor)
	{
		
		//initial the variables for air temperature
		float Hn  <- list_sunrise at (nb_days); //daily sunrise time
		float Ho <- list_sunset at (nb_days); //daily sunset time
		float Hp <- list_sunrise at (nb_days+1); //next day sunrise time
		float Tn  <-  list_daily_min_temp at (nb_days);
		float Tx  <- list_daily_max_temp at (nb_days);
		float Tp <- list_daily_min_temp at (nb_days + 1);	
		float Hx <- 12.0; // set the maximum temperature of each day to noon (12 in pacific standard time)
		float To <- Tx - 0.39 * (Tx - Tp);
		float alpha  <- Tx - Tn;
		float R  <- Tx - To;
		float b <- (Tp - To)/ sqrt(abs(Hp - Ho));
		float pi <- 3.14;
		float HourAvgTempC <- 0.0;
		prevAirTempHour <- fAirHourAvgTempK - 243.15;


		if (calTime > Hn and time <= Hx) { HourAvgTempC <-Tn + (alpha * (((calTime - Hn)/ (Hx - Hn)) * (pi/2)));}
		else if (calTime > Hx and time < Ho)  { HourAvgTempC <- To + (R* sin(((pi/2) + (((calTime - Hx)/4) * (pi/2)))));}
		else if (calTime > Ho or time <= Hp) { HourAvgTempC <- To + (b * (sqrt(abs(calTime- Ho))));}

		return HourAvgTempC + 273.15;	//return the temperature in degre Kelvin
	}
	
	
	
	/**************************************************************
	 * compute for everyday the statistics of all Aedes activities
	 **************************************************************/
	reflex compute_duration_stat when: bCheckActStat and cycle > 1 and  bNewDay
	{
		float t <- bVarComputeExecTime ? machine_time : 0; 
		int nb <- 0;
		iNbAvgDailyTb <- length(Aedes) > 0? int(sum(Aedes collect (each.iCmpDailyTb))/length(Aedes)): 0;
		iNbAvgDailyTn <- length(Aedes) > 0? int(sum(Aedes collect (each.iCmpDailyTn))/length(Aedes)): 0;
		iNbAvgDailyBs<- length(Aedes) > 0? int(sum(Aedes collect (each.iCmpDailyBs))/length(Aedes)): 0;
		iNbAvgDailyShd<- length(Aedes) > 0? int(sum(Aedes collect (each.iCmpDailyShd))/length(Aedes)): 0;
		iNbAvgDailyFt<- length(Aedes) > 0? int(sum(Aedes collect (each.iCmpDailyFt))/length(Aedes)): 0;
		iNbAvgDailyWd <- length(Aedes) > 0? int(sum(Aedes collect (each.iCmpDailyWd))/length(Aedes)): 0;
		fAvgDailyFlyDist <- length(Aedes) > 0? float(sum(Aedes collect (each.fDailyFlightDist))/length(Aedes)): 0;	
		
		list<float> list_blood <- Aedes collect ("blood" in each.map_actTotDailyDur.keys ? each.map_actTotDailyDur["blood"] : 0.0);
		list<float> list_layEgg <- Aedes collect ("bs" in each.map_actTotDailyDur.keys ? each.map_actTotDailyDur["bs"] : 0.0);
		list<float> list_shade <- Aedes collect ("shade" in each.map_actTotDailyDur.keys ? each.map_actTotDailyDur["shade"] : 0.0);
		list<float> list_nectar <- Aedes collect ("nectar" in each.map_actTotDailyDur.keys ? each.map_actTotDailyDur["nectar"] : 0.0);
		list<float> list_flyTo <- Aedes collect ("flyTo" in each.map_actTotDailyDur.keys ? each.map_actTotDailyDur["flyTo"] : 0.0);
		list<float> list_wander <- Aedes collect ("wander" in each.map_actTotDailyDur.keys ? each.map_actTotDailyDur["wander"] : 0.0);
		save[nb_days, list_blood] type: csv to: "blood-duration.csv";

		//duration for a blood meal
		list_blood <- list_blood sort_by (each);

		stat_takeBlood[0] <- first(list_blood);
		stat_takeBlood[1] <- last(list_blood);
		stat_takeBlood[2] <- mean(list_blood);
		nb <- length(list_blood);
		stat_takeBlood[3] <- list_blood[int(nb/4)];
		stat_takeBlood[4] <- list_blood[int(nb/2)];
		stat_takeBlood[5] <- list_blood[int(3*nb/4)];
		
		
		// duration for a resting
		list_shade<- list_shade sort_by (each);
		stat_rest[0] <- first(list_shade);
		stat_rest[1] <- last(list_shade);
		stat_rest[2] <- mean(list_shade);
		nb <- length(list_shade);
		stat_rest[3] <- list_shade[int(nb/4)];
		stat_rest[4] <- list_shade[int(nb/2)];
		stat_rest[5] <- list_shade[int(3*nb/4)];
		 
		
		//
		ask Aedes {
			loop act over: list_KeyTargetName + ["flyTo", "wander"]{//["flyTo", "wander", "waitEggMat"] {
				map_actTotDailyDur[act] <- 0.0;
			}
			
		// Re-init counter
		fDailyFlightDist <- 0.0;//save it to print		
		iCmpDailyTb <- 0;
		iCmpDailyTn <- 0;
		iCmpDailyBs <- 0;
		iCmpDailyShd <- 0;
		iCmpDailyFt <- 0;
		iCmpDailyWd <- 0;								
		}	
		
		if (bVarComputeExecTime) {t1 <- t1 + machine_time - t;	}
			
	}
	
	
	/*********************************************************
	 * World update all Aedes stocks here for display chart
	 *********************************************************/
	reflex update_stock when: bUpdateStck and (cycle > 1 or bNewDay)
	{
		float t <- bVarComputeExecTime ? machine_time : 0;  
		totEggs <- sum(SpatObj collect (each.fTotStckE));
		totLarvae <- sum(SpatObj collect (each.fTotStckL));
		totPupae <- sum(SpatObj collect (each.fTotStckP));
		totNewAdults <- sum(SpatObj collect (each.iNbNewFemaleA));
		if (bVarComputeExecTime) {t2 <- t2 + machine_time - t;	}	
	}
	

}


/////////////////////////////////////////////////////////
// 				ALL OF SPECIES
/////////////////////////////////////////////////////////
entities 
{
	
	/*****************************************************************************************************************************
	 * SPECIE DISPERSAL 
	 * (use to create a shapefile represent a BS location and the maximum dispersal of Aedes population)
	 ****************************************************************************************************************************/	
	species Dispersal
	{
		float fMaxRadius;
		string strOrigSpatObj;
		rgb color;
		
		init
		{
			fMaxRadius <- 0.0;
			color <-rgb("blue");
			ask one_of(Aedes) 
			{
				myself.location <- bornLoc;
			}
		}
		
		aspect default
		{
			draw circle(fMaxRadius) color: color;
		}
	}
	
	
	
	/*******************************************************
	 * 			SPECIE AEDES
	 ******************************************************/
	species Aedes skills:[moving]
	{
		string strState;	
		bool bActive;
		bool bGonoBegan;
		float fCurrGonoPhysioDev;
		int iAge; //is initiated by the SpatObj when new Agent Aedes has been created
		bool bCurrActDone;
		float fTimeStartResting;
		bool bBestTargetReached;
		string strTarget;
		int iCmpTimeflyTo;
		float fTimeNeedAct;
		float fMinBloodNeedDevEgg;
		int iCmpGonoCycle;
		SpatObj bsToOvipo;
		bool bMated;
		SpatObj bestTarget;
		point bestTargetPoint;		
		int iNbStates; //private
		int iNbTargets; // private
		bool bMovingTrack;
		list<point> movingPts_list;
		geometry trail;
		
		float fCurrEnergy;
		float fQntStckBlood;
		float fEnergyGainInAct <- 0.0;
		float fEnergyLostInAct <- 0.0;		
		list<float> list_finalTargertNeedWg <- [0.0,0.0,0.0,0.0,0.0]; // final need weight for each target
		
		geometry perception_area;// update: cone({heading-30,heading+30}) inter circle(PERCEPTION_RADIUS);
		SpatObj mySpatObj;// update: SpatObj closest_to self;//first(SpatObj overlapping(self));
		list<SpatObj> list_limit <- [] ;
		
		bool bNeedBloodGono;
		bool bPhysioDevDone;
		float fTimeLeftLayStck <- 0.0;
		
		list<list> dailyAct <- [];
		map<string,float> map_actTotDailyDur; //save total duration for each action (6 actions) for an Aedes
		
		list act_leftPrevDay <- []; //temporary var
		list last_actDay <- [];
		
		int iDayNeedActif;
		
		////////////// statistic variables ///////////////////////
		SpatObj bornObj;	
		point bornLoc; //Aedes born location
		float fMaxLifeDisps; // maximum Dispersal distance
		float fCoordX;
		float fCoordY;
		
		float fDailyFlightDist;
		int iCmpDailyTb;
		int iCmpDailyTn;
		int iCmpDailyBs;
		int iCmpDailyShd;
		int iCmpDailyFt;
		int iCmpDailyWd;
		int iCmpAct; // count the number of activities in the hold life of Aedes
		int iCmpDailyTrail; // number of flight per day
		
		//----------------------------------
		// init
		//----------------------------------
		init
		{
			speed <- (rnd((MAX_FLIGHT_SPEED - MIN_FLIGHT_SPEED)*1000) / 1000 + MIN_FLIGHT_SPEED); //GAMA pseudo-variable (default m/s) = 0.5 - 1 m/s or 30 - 60m/mn !!!!!!!! CHECK if work automatically?		
			iDayNeedActif <- int(rnd(3*1000)/1000);//betwenne 0 - 2 days resting after emerging
			fMinBloodNeedDevEgg <- gauss(AVG_MIN_BLOOD_NEED_DEVEGG, ((UPPERLIMIT_MIN_BLOOD_NEED_DEVEGG - LOWERLIMIT_MIN_BLOOD_NEED_DEVEGG)/6));
			iNbStates <- length(mapStateIndex);
			iNbTargets <- length(mapTargetIndex);
			bActive <- true;
			bGonoBegan <- false;
			fCurrGonoPhysioDev <- 0.0; 
			fDailyFlightDist <- 0.0;
			bCurrActDone <- true;
			bBestTargetReached <- false;
			bestTarget <- nil;
			strTarget <- "";
			fTimeNeedAct <- 0.0;
			fQntStckBlood <- 0.0;
			iCmpGonoCycle  <- 0;
			bsToOvipo <- nil;
			bMated <- flip(MATING_PROB); 
			bestTarget <- nil;
			bMovingTrack <- false;
			movingPts_list <- [];
			bNeedBloodGono <- false;
			fTimeLeftLayStck <- 0.0;
			fMaxLifeDisps <- 0.0;

			bPhysioDevDone <- false;
			fCurrEnergy <- MAX_ENERGY;			
			//initial state
			if(!bMated) 
			{
				strState <- "virgin";
			} 
			else if (fQntStckBlood <= 0.0) //bMated and never feed blood
			{
				strState <- "ovipo";
				fCurrGonoPhysioDev <- 0.0;
			} 
			else 
			{
				strState <- "gono";
			}

			do updateCurrEnergy;
			do updateTargetWeight;

		}
		

		//---------------------------
		// Reflexes
		//---------------------------

		/**********************************************************
		 * R1: updateGonoPhysioDev
		 * update for each hour the gonotrophic development rate
		 * and decrease blood stock (assume to use blood for dev. eggs)
		 * "hourlyPhysioDevGono" fluctuate between [0, 1]
		 /
		 ***********************************************************/
		reflex updateGonoPhysioDev when: bGonoBegan and bNewHour and !bPhysioDevDone
		{
			float t <- bVarComputeExecTime ? machine_time : 0; 
			if(bNeedBloodGono = false) // check if there is not enough blood to dev eggs
			{
				hourlyPhysioDevGono <- 0.0;
				if(avgHourAirTempK_matrix[iSimRealDay, int(curTime)] != 0)//in case of divide by 0
				{
					hourlyPhysioDevGono <- RO25_GONO* ((avgHourAirTempK_matrix[iSimRealDay, int(curTime)]/298)*exp(DHA_GONO/R*(1/298 - 1/avgHourAirTempK_matrix[iSimRealDay, int(curTime)])))/ (1 + exp(DHH_GONO/R*(1/THALF_GONO - 1/avgHourAirTempK_matrix[iSimRealDay, int(curTime)])));
				}
				fCurrGonoPhysioDev <- fCurrGonoPhysioDev + hourlyPhysioDevGono;
				fQntStckBlood <- fQntStckBlood - BLOOD_GONOUSE_ULPERHOUR;//ul - eliminate blood stock in order to use for egg dev (cf: Christopher, 1960 => in general 3mg (full meal) of blood is eliminate withing 10h)
				
				if(fQntStckBlood <= fMinBloodNeedDevEgg)
				{
					bNeedBloodGono <- true;
				}					
			}
			
			do checkPhysioDevDone;	
			if (bVarComputeExecTime) {t4 <- t4 + machine_time - t;	}
		}

		/**********************************************************
		 * R2: updateControlParam
		 * update all the control variables for each day (midnight):
		 * - age is incresed each day 
		 * - daily flight distance reinitiate to zero
		 **********************************************************/
		reflex updateControlParam when: bNewDay
		{
			float t <- bVarComputeExecTime ? machine_time : 0; 		
			iAge <- iAge + 1;
			if iAge >= MAX_AGE_ADULT 
			{
				do killMe("age");
			}
			if (bVarComputeExecTime) {t14 <- t14 + machine_time - t;	}	
			 			
		}		
		
		/**********************************************************
		 * R3: checkCurrentAct
		 *********************************************************/
		reflex checkCurrentAct when: bActive
		{
			float t <- bVarComputeExecTime ? machine_time : 0; 
			
			mySpatObj <- first(SpatObj overlapping self);
			if (mySpatObj = nil) {mySpatObj <- SpatObj closest_to self;}
			
			if(bShowLimitArea)
			{
				ask(SpatObj) 
				{
					ask (world) {myself.color <-  adjustObjColor(myself.spcClass.color);}
				}	//reinitiate the color
			}
			
			//---------------------------------------------
			if(bCurrActDone)//previous activity done
			{
				if(iAge < iDayNeedActif) // if a female is a new born, stay rest for 1-2 days!!
				{
					do rest;
				}
				else
				{
					do decideBestAct;	
				}
			}
			else //continue the activity
			{
				switch strTarget
				{
					match "nectar" {do takeNectar;}
					match "blood" {do takeBlood;}
					match "bs" 
					{
						do layEgg; //check welcomed capacity will do in action layEgg
					}
					match "shade"
					{
						do rest;
					}
					default //case of previous target = nothing
					{
						point prv_loc <- location;
						//do update_limit;
						do wander speed: speed/step bounds: mySpatObj.geom_limit;
						if(bCheckActStat)
						{
							do save_activity("wander",step);
						}
						do saveFlightDistance(prv_loc); //to memorise max flight
						strTarget <- "";
						bestTarget <- nil;
						bBestTargetReached <- false;
					}
				}
			}
			
			
			if(bCheckFinalFuzzyScore)
			{	//save a file for each aedes
				save [name, nb_days, curTime, list_finalTargertNeedWg[0], list_finalTargertNeedWg[1], list_finalTargertNeedWg[2], list_finalTargertNeedWg[3], list_finalTargertNeedWg[4], strTarget] type: "csv" to: "fSTarget-"+name+".csv";
			}
			if (bVarComputeExecTime) {t5 <- t5 + machine_time - t;}	
		}
		

		
		/*******************************************************************
		 * Kill me: kill Aedes agent and save all needed activity to a file
		 *******************************************************************/
		action killMe(string cause)
		{
			if(bTrackAedes)
			{
				trail <- line(movingPts_list);
				create dailyTrail  with: [aedesName::name, shape::trail];
				release movingPts_list;	
			}
			if(bCheckMaxDisp)
			{
				string tmp <- string(bornObj) +","+ name + ","+ iAge+ ","+ fMaxLifeDisps+ ","+ fCoordX+ "," + fCoordY + ","+ cause + "\n";
				txtSaveDeath <-txtSaveDeath + tmp;
				// save the maximum location of an Aedes
				//create dispersal_loc with: [aedesName::name, sim_day::nb_days, shape::point(fCoordX, fCoordY), coorX::fCoordX, coorY::fCoordY];
				//save dispersal_loc type:"shp" to: dispersShpName with:[aedesName::"AEDES", sim_day::"DAY_SIM", coorX::"coorX", coorY :: "coorY"];
			}
			if ( bCheckMaxDisp){ask bornObj { do updateDispersalCircle(myself);}}
			do die;
		}
		
		/*******************
		 * AB: decideBestAct
		 * is called when an activity finished from the previous step
		 * new female must stay rest 1- 2 days before being active
		 *
		 *******************/
		action decideBestAct
		{
		
			if(strState = "virgin")
			{
				if(flip(MATING_PROB))
				{
					strState <- "ovipo";
					do updateCurrEnergy;
					do updateTargetWeight; //shade 0.5 blood 1 nectar 0.2 bs 0
				}
			}
				
			
			if(bestTarget = nil) //au commencement et Ã  la fin d'une activitÃ©
			{
					if(fuzzyChsBestTarget()) 
					{
						do flyToBestTarget;
					}
			}
			else //if bestTarget has a value (!= nil)
			{			 	
				if(bBestTargetReached)
				{
					if(applyAttempSurvivalRate())
					{
			
							switch strTarget
							{
								match "nectar" 
								{
									fTimeNeedAct <- getTimeNeed(); 
									time_nectar <- fTimeNeedAct;
									
									if(bCheckActStat)
									{
										do save_activity(strTarget,time_nectar);
									} 
									fEnergyGainInAct <- 0.0;
									fEnergyLostInAct <- 0.0;
									do takeNectar;
								}
								match "shade" 
								{
									fEnergyGainInAct <- 0.0;
									fEnergyLostInAct <- 0.0;
									do rest; //average 1- 3 days in second
								}
								match "bs" 
								{								
									if (bPhysioDevDone and checkEnvirToLayEgg() ) 
									{
										bsToOvipo <- bestTarget;
										if(fTimeLeftLayStck > 0.0)
										{
											fTimeNeedAct <- fTimeLeftLayStck;
										}
										else
										{
							                fTimeNeedAct <- getTimeNeed(); 
										}
										
										time_lay_egg <- fTimeNeedAct;
										if(bCheckActStat)
										{
											do save_activity(strTarget,time_lay_egg);
										} 
										fEnergyGainInAct <- 0.0;
										fEnergyLostInAct <- 0.0;
										
										ask mySpatObj
										{
										 	iHourLayEgg	<- iHourLayEgg +1;
										}
										do layEgg;
									}
									else //libÃ©rer le bestTarget pour le nouveau dÃ©cision
									{
										point prv_loc <- location;
										if(debug){write "wander for bs";}
										//do update_limit;
										do wander speed: speed/step bounds:  mySpatObj.my_neighbors; //geom_limit; /*bounds: gDailyFlightCircle;*/
										if(bCheckActStat)
										{
											do save_activity("wander",step);
										} 	
										do saveFlightDistance(prv_loc);				
										strTarget <- "";
										bestTarget <- nil;//in order to lauch a new fuzzy decision
										bBestTargetReached <- false; //use after new fuzzy decision
										fEnergyGainInAct <- 0.0;
										fEnergyLostInAct <- step * WANDER_ENERGLOST_RATE;
										do updateCurrEnergy;
										do updateTargetWeight;	
									}
								}
								match "blood" 
								{
									if (flip(BITING_SUCC_RATE)) 
									{
										fTimeNeedAct <- getTimeNeed(); 
										time_stat_takeBlood <- fTimeNeedAct;									
										if(bCheckActStat)
										{
											do save_activity(strTarget,time_stat_takeBlood);
										} 
										fEnergyGainInAct <- 0.0; //initiate the value before the action
										fEnergyLostInAct <- 0.0;
										if(bCheckNbBiteOvipoPerHourInObj)
										{
											ask mySpatObj
											{
											 	iHourBite	<- iHourBite +1;
											}
										}
										do takeBlood;
									} 
									else 
									{
										point prv_loc <- location;
										//do update_limit;
										do wander speed: speed/step bounds: mySpatObj.geom_limit;  /*bounds: gDailyFlightCircle;*/			
										if(bCheckActStat)
										{
											do save_activity("wander",step); 
										}
										do saveFlightDistance(prv_loc);
																				
										strTarget <- "";
										bestTarget <- nil;//lauch a new fuzzy decision
										bBestTargetReached <- false; //use after new fuzzy decision
			
										fEnergyGainInAct <- 0.0;
										fEnergyLostInAct <- step * WANDER_ENERGLOST_RATE;
										do updateCurrEnergy;
										do updateTargetWeight;	
									}
								}
								default //include strTarget = "not"
								{
										point prv_loc <- location;
										//do update_limit;
										do wander speed: speed/step bounds: mySpatObj.geom_limit;  /*bounds: gDailyFlightCircle;*/
										if(bCheckActStat)
										{
											do save_activity("wander",step); 
										}
										do saveFlightDistance(prv_loc);
																				
										strTarget <- "";
										bestTarget <- nil;//lauch a new fuzzy decision
										bBestTargetReached <- false; //use after new fuzzy decision
										fEnergyGainInAct <- 0.0;
										fEnergyLostInAct <- step * WANDER_ENERGLOST_RATE;
										do updateCurrEnergy;
										do updateTargetWeight;
								}									
							}//end switch				
					} //end if attempSurv!
					else //die while attemping to bite
					{
						do killMe("bite");
					}
				} //end bestTargetReached
				else //if bestTarget != nil and !bBestTargetReached
				{
					do flyToBestTarget; //continue to find a target
				}
			}
		}
		
		
		/*************************************************************************************
		 * AA: ask the spatial object where the mosquito is the capabability having new eggs
		 ************************************************************************************/
		bool checkEnvirToLayEgg
		{
			ask (mySpatObj)
			{
				if(getCapaHaveEgg() <= 0.0)
				{
					return false; // Aedes cannot lay egg on it selected BS
				}
				else
				{					
					return true; //Aedes can lay egg inside
				}
			}	
		}

		/*******************
		 * A1: rest 
		 * this action is called by an agent SpatObj when the light is considering as non-favorable
		 *
		 *******************/
		action rest
		{
			fTimeStartResting  <- time;	
			bActive <- false;
		}

		/*******************
		 * A2: wake up and then fly randomly for this step before deciding a new act
		 *******************/		
		action wakeUp
		{
			if(strState != "gono") or (strState ="gono" and bNeedBloodGono) or (strState = "gono" and bPhysioDevDone) 
			{
				if(bCheckActStat)
				{
					do save_activity(strTarget, time - fTimeStartResting );
				}
				//release or the conditions to lauch a new decision making
				fTimeStartResting <- 0.0;			
				bBestTargetReached <- false;
				bestTarget <- nil;
				strTarget <- "";
				bActive  <- true;	
				bCurrActDone <- true; //finish resting
				
				// ------ fly randomly after waking up ------
				point prv_loc <- location;
				//do update_limit;
				do wander speed: speed/step bounds: mySpatObj.geom_limit;	
				if(bCheckActStat)
				{
					do save_activity("wander",step); 
				}
				do saveFlightDistance(prv_loc);
				do updateCurrEnergy;
				do updateTargetWeight; //normally, the mosquito must want less resting desire because it earns energy while resting				
			}
		}
		
		
		/*******************
		 * A3: take blood 
		 * return fQntStckBlood: use as a base element to calculate current energy level 
		 *******************/
		action takeBlood
		{
			float fQntBloodFeed <- 0.0;
			
			if(fTimeNeedAct <= step)//activity is finished in a step
			{
				bCurrActDone <- true;
				strTarget <- "";
				bBestTargetReached <- false;
				bestTarget <- nil;
				bNeedBloodGono <- false; //in case of need blood during gonotrophic cycle
				fQntBloodFeed <- (rnd(BITING_SPEED_UL_SEC*10000)/10000) * fTimeNeedAct;
				fQntStckBlood <- fQntStckBlood + fQntBloodFeed; //need to calculate for each step because this variable has an impact to the beginning of Gonotrophique cycle
				fEnergyGainInAct <- fEnergyGainInAct + (BLOOD_GAIN_ENERGY_SPEED*fTimeNeedAct);
				fEnergyLostInAct <- fEnergyLostInAct + (fTimeNeedAct * TAKEBLOOD_ENERGLOST_RATE);				
				do updateCurrEnergy;
				do updateTargetWeight;
			}
			else
			{
				bCurrActDone <- false;
				fTimeNeedAct <- fTimeNeedAct - step; //calculated time left to do an activity
				fQntBloodFeed <- (rnd(BITING_SPEED_UL_SEC*10000)/10000) * step;
				fQntStckBlood <- fQntStckBlood + fQntBloodFeed; //need to calculate for each step because this variable has an impact to the beginning of Gonotrophique cycle	
				fEnergyGainInAct <- fEnergyGainInAct + BLOOD_GAIN_ENERGY_SPEED*step;
				fEnergyLostInAct <- fEnergyLostInAct + step * TAKEBLOOD_ENERGLOST_RATE;
			}
			
			//check for the beginning of gonotrophic cycle only for "ovipo" female (it depends on the level of blood fed)	
			//the virgin female doen't concern with this process cause the blood meal won't change anything to her.
			if(strState = "ovipo" and fQntStckBlood >= fMinBloodNeedDevEgg )
			{
				bGonoBegan <- true; 
				strState <- "gono";
				iCmpGonoCycle<- iCmpGonoCycle + 1;
			}			
		}
		
		/*******************
		 * A4: take nectar 
		 *******************/		
		action takeNectar
		{
			if(fTimeNeedAct <= step)
			{
				bCurrActDone <- true;
				strTarget <- "";
				bBestTargetReached <- false;
				bestTarget <- nil;
				fEnergyLostInAct <- fEnergyLostInAct + TAKENECTAR_ENERGLOST_RATE* fTimeNeedAct;
				fEnergyGainInAct <- fEnergyGainInAct + NECTAR_GAIN_ENERGY_SPEED*fTimeNeedAct;
				do updateCurrEnergy; //mn
				do updateTargetWeight;
			}
			else
			{
				bCurrActDone <- false; //indeed, don't need this process cause the variable has already assigned a value false since finishd flyTotarget
				fTimeNeedAct <- fTimeNeedAct - step; //calculated time left to do an activity
				fEnergyLostInAct <- fEnergyLostInAct + TAKENECTAR_ENERGLOST_RATE* step;
				fEnergyGainInAct <- fEnergyGainInAct + TAKENECTAR_ENERGLOST_RATE* step;
			}		
		}
		
		/*******************
		 * A5: lay eggs
		 *******************/			
		action layEgg
		{
			int eggsLaid <- 0;
			float fBloodUsedLayEgg <- 0.0;
			float fSpeed <- (rnd((MAX_SPEED_LAYEGG - MIN_SPEED_LAYEGG)*10000)/10000 + MIN_SPEED_LAYEGG);
		 	bsToOvipo <- mySpatObj;

			if(fTimeNeedAct <= step)//activity is finished in a step
			{
				bCurrActDone <- true;
				strTarget <- "";
				strState <- "ovipo";
				bBestTargetReached <- false;
				bestTarget <- nil;
			 	eggsLaid <- int(fSpeed * fTimeNeedAct); //calculate the nb of eggs to lay for a step (in sec). Aedes can lay 1-3 eggs/mn
				fBloodUsedLayEgg <- eggsLaid * BLOOD_USE_LAYEGG_RATE; // eggs * (mg/eggs) = mg
				fQntStckBlood <- fQntStckBlood - fBloodUsedLayEgg; //decrease stock of blood that is used in laying eggs.
				fEnergyLostInAct <-  ENERGY_LOST_LAYEGGE_SEC * fTimeNeedAct;

				do updateCurrEnergy; //mn
				do updateTargetWeight;//(increase blood desire weight, decreas BS desire one)
			}
			else 
			{
				bCurrActDone <- false;
				fTimeNeedAct <- fTimeNeedAct - step; //calculated time left to do an activity
			 	eggsLaid <- int(fSpeed * step); //calculate the nb of eggs to lay for a step (in sec). Aedes can lay 1-3 eggs/mn 
				fBloodUsedLayEgg <- eggsLaid * BLOOD_USE_LAYEGG_RATE; // eggs * (mg/eggs) = mg
				fQntStckBlood <- fQntStckBlood - fBloodUsedLayEgg; //decrease stock of blood that is used in laying eggs.
				fEnergyLostInAct <- ENERGY_LOST_LAYEGGE_SEC * step;
			}
			

			//------ crate new Aedes agent
			 if(bsToOvipo != nil)
			 {
				
				// add new eggs to the stock of eggs in the container where oviposition is done.			
				ask bsToOvipo
				{
					float fBsInRate <- 0.5;
					bool bLayInside <- true;
					//calculate the rate to lay eggs inside based on the current proportion of water-filled breeding site in on the total (in and out) of the SpatObj
					if(iNbBsInWithWater + iNbBsOutWithWater > 0)
					{
						fBsInRate <- iNbBsInWithWater/(iNbBsInWithWater + iNbBsOutWithWater);
					}
				 	
				 	//-- choose where to lay egg (in or out) based on the bsInRate
				 	if(!flip(fBsInRate)) 
				 	{
						bLayInside <- false; //lay outside
					}

					//update the stock of eggs waiting for water
					if(eggsLaid < getCapaHaveEgg())//Aedes can lay all eggs
					{
					 	if(bLayInside)
					 	{
					 		iWaitStckEggIn  <- iWaitStckEggIn + eggsLaid;
					 		if(bCheckStck and eggsLaid != 0)
					 		{
								string tmp <- ""+ nb_days + ","+ myself.name + ","+ self.name+ ","+ self.strClass + "," + self.fSpatObjWaterTempCIn+ ","+ eggsLaid+ ",in"+ "\n";
								txtEggLaid <-txtEggLaid + tmp;
							}		 			
					 	}
					 	else	
					 	{
					 		iWaitStckEggOut  <- iWaitStckEggOut + eggsLaid;
					 		if(bCheckStck and eggsLaid != 0)
					 		{
					 			string tmp <- ""+ nb_days + ","+ myself.name + ","+ self.name+ ","+ self.strClass + "," + self.fSpatObjWaterTempCOut+ ","+ eggsLaid+ ",out"+ "\n";
								txtEggLaid <-txtEggLaid + tmp;
					 		}		 			
					 	}
					 	myself.fTimeLeftLayStck <- 0.0;		
					}
					else // lay only in the limit of the capacity having eggs
					{
						if(bLayInside)//rate to choose bs inside
				 		{
							iWaitStckEggIn <-  iWaitStckEggIn + iCapaHaveEgg;
					 		if(bCheckStck and eggsLaid != 0)
					 		{
					 			string tmp <- ""+ nb_days + ","+ myself.name + ","+ self.name+ ","+ self.strClass + "," + self.fSpatObjWaterTempCIn+ ","+ iCapaHaveEgg + ",in"+ "\n";
								txtEggLaid <-txtEggLaid + tmp;
					 		}		 			
						}
						else //lay egg to BS OUTSIDE
				 		{
					 		iWaitStckEggOut <-  iWaitStckEggOut + iCapaHaveEgg;
					 		if(bCheckStck and eggsLaid != 0)
					 		{
					 			string tmp <- ""+ nb_days + ","+ myself.name + ","+ self.name+ ","+ self.strClass + "," + self.fSpatObjWaterTempCOut+ ","+ iCapaHaveEgg + ",out"+ "\n";
								txtEggLaid <-txtEggLaid + tmp;
					 		}
				 		}
					 		
					 	myself.fTimeLeftLayStck <- myself.fTimeNeedAct + ((eggsLaid - iCapaHaveEgg)/fSpeed);
					 	myself.list_finalTargertNeedWg[mapTargetIndex["bs"]] <- 1.0; //force to find a new breeding site
					 	myself.bCurrActDone <- true;
						myself.strTarget <- "";
						myself.strState <- "ovipo";
						myself.bBestTargetReached <- false;
						myself.bestTarget <- nil;
						ask myself {do updateCurrEnergy;} //mn
						ask myself {do updateTargetWeight;}//(increase blood desire weight, decreas BS desire one)
					 }	
				}
			}			
		}
	
		/**************************************************************************
		 * choose the best target based on need and availibility of target in a SpatObj 
		 * ************************************************************************/
		bool fuzzyChsBestTarget
		{	
			float t <- bVarComputeExecTime ? machine_time : 0;  
			list<SpatObj> list_detectedObjs <- [];
			float fDistToObj <- 0.0;
			float fAttValue <- 0.0;
			float fTotScore <- 0.0;
			bool bTargetFound <- false;
			float fMaxDistMyObj <- 1.0;
			list<float> listSelDist <- [];
			list<point> listSelPoint <- [];
			
			////// create a perception_cone
			point p1 <- {location.x + PERCEPTION_RADIUS * cos(heading+30),location.y + PERCEPTION_RADIUS * sin(heading+30) };
			point p2 <- {location.x + PERCEPTION_RADIUS * cos(heading-30),location.y + PERCEPTION_RADIUS * sin(heading-30) };
			perception_area <- polygon([location,p1,p2]);
			
			/////// detect the spatial objects inside the perception area /////////
			list_detectedObjs <- mySpatObj.my_neighbors where (each overlaps perception_area);
			if(mySpatObj  != nil) 
			{
				add mySpatObj to: list_detectedObjs;
			}
			
			
			////////////////// calculate the final score of each target on an spatial obj /////////
			if(!empty(list_detectedObjs))
			{
				float fTotDensDecHuman <- 0.0; //use to count the total number of human pop in perceived objets
				matrix<float> mat_finalAttrctScore <- 0.00 as_matrix({length(list_detectedObjs),iNbTargets}); // attraction score of each target regarding to each object spatial		
				loop i from: 0 to: length(list_detectedObjs) - 1 // for all detected objs
				{
					
					point pt <- location; //initialisation// point ((list_detectedObjs[i] closest_points_with location) at 0); 
					//point pt <- any_location_in(list_detectedObjs[i]); // choose randomly one location in a selected obj
					geometry interZone <- (list_detectedObjs[i].shape inter perception_area);
					
					//if Aedes reaches the border of simulation zone, turn 90° till finding something
					loop while: interZone = nil 
					{
					    heading <- heading - 90;
					    perception_area <- cone({heading-30,heading+30}) inter square(PERCEPTION_RADIUS);
						interZone <- (list_detectedObjs[i].shape inter perception_area);
					}
					pt <- any_location_in(interZone);			
					fDistToObj <- self distance_to pt;
							
					fMaxDistMyObj <- self distance_to point ((mySpatObj closest_points_with location) at 0);

					if(fDistToObj < DIST_NOEFF) 
					{
						fDistToObj <- DIST_NOEFF; //IF THE OBJ IS IN THE LIMIT OF MY OBJ, NO EFFECT OF DIST TO THE DECISION (CAUSE DIVIDE BY 1)
					}
				    else if(fDistToObj < fMaxDistMyObj) 
					{
						fDistToObj <- fMaxDistMyObj; //IF THE OBJ IS IN THE LIMIT OF MY OBJ, NO EFFECT OF DIST TO THE DECISION (CAUSE DIVIDE BY 1)
					}	
				
					add fDistToObj to: listSelDist;
					add pt to:listSelPoint;
					
					loop j from: 0 to: iNbTargets -1 // for all targets
					{
						
						ask list_detectedObjs[i] 
						{
							fAttValue <- list_attrctTarget[j]; 
							if(list_KeyTargetName[j] = "blood") //for blood, fAttValue represente a density of human in an obj
							{
								fTotDensDecHuman <- fTotDensDecHuman + fAttValue; // accumulate the density of human from each detected obj
							}
						}
						mat_finalAttrctScore[i,j] <- (list_finalTargertNeedWg[j] * fAttValue)/ (fDistToObj^0.1); //final score for each obj and each target
						
						//make a total score for all the target in each detected obj except for blood target
						if(list_KeyTargetName[j] != "blood") 
						{
							fTotScore <-  fTotScore + mat_finalAttrctScore[i,j];
						}
					}	
					
				}
				
				////// Reassigne the blood attraction value to perceiving obj ////////////
				if(fTotDensDecHuman > 0.000)
				{
					loop i from: 0 to: length(list_detectedObjs) - 1 // for all detected objs
					{
						int j <-mapTargetIndex["blood"];
						ask list_detectedObjs[i] 
						{
							fAttValue <- list_attrctTarget[j]; //= density of human in an obj i
						}
						
						mat_finalAttrctScore[i,j] <- (list_finalTargertNeedWg[j] * (fAttValue/fTotDensDecHuman))/ (fDistToObj^0.1);
						fTotScore <-  fTotScore + mat_finalAttrctScore[i,j];
						
						// ---------- for each detected obj i, save the selection process if ask for it!----------------
						if(bTrackFuzzyProcess)
						{
							save["nectar", list_detectedObjs[i].name, fDistToObj, list_finalTargertNeedWg[mapTargetIndex["nectar"]], mat_finalAttrctScore[i,mapTargetIndex["nectar"]]] type: csv to: fileFuzzyProc;
							save["blood", list_detectedObjs[i].name, fDistToObj, list_finalTargertNeedWg[mapTargetIndex["blood"]], mat_finalAttrctScore[i,mapTargetIndex["blood"]]] type: csv to: fileFuzzyProc;
							save["shade", list_detectedObjs[i].name, fDistToObj, list_finalTargertNeedWg[mapTargetIndex["shade"]], mat_finalAttrctScore[i,mapTargetIndex["shade"]]] type: csv to: fileFuzzyProc;
							save["bs", list_detectedObjs[i].name, fDistToObj, list_finalTargertNeedWg[mapTargetIndex["bs"]], mat_finalAttrctScore[i,mapTargetIndex["bs"]]] type: csv to: fileFuzzyProc;
						}
					}
				}

				
				/////////// choose the best target with stochastic choise on the final attract score (mat_finalAttrctScore) //////////
				float fChooseScore <- rnd(fTotScore * 1000)/1000;
					
				int i <- 0;
				int j <- 0;
				float fTmpScore <- 0.0;//mat_finl[0,0]
				//float fTmpScore <- mat_finalAttrctScore[i,j];//mat_finl[0,0]
							
				loop while: i <= (length(list_detectedObjs) - 1)	
				{
					loop while: j <= (iNbTargets-1)
					{
						fTmpScore <- fTmpScore + mat_finalAttrctScore[i,j];	
						if(fChooseScore <= fTmpScore)
						{
							bestTarget <- list_detectedObjs[i];
							//verify if the target can be reached because of the porosity of the actual object
							bool bProbPass <- flip((mySpatObj.fPorosity) * (bestTarget.fPorosity));
							if(bProbPass) // if Aedes can g through the actual obj and reach the target, assign the point of the bestTarget
							{
								bestTargetPoint <- listSelPoint[i];
							}
							else //assign the border of mySpatObj and Aedes can't find its best target (so strTarget = "")
							{
								bestTargetPoint <- point ((mySpatObj closest_points_with location) at 0);
								//bestTargetPoint <- mySpatObj farthest_point_to(self.location);
								strTarget <- ""; 
							}

							strTarget  <- list_KeyTargetName[j];
							if(bTrackFuzzyProcess)
							{
								float estimFlyDist <-  self distance_to bestTargetPoint;
								save [nb_days, hour_of_a_day, name, mySpatObj, bestTarget, strTarget,estimFlyDist, listSelDist[i]] type: csv to: fileSelBestTarget;
							}	
							
							i <- 999; //to stop the loop
							j <- 999; // to stop the loop
																	
						}
						j <- j + 1;
					}
					i <- i + 1;
					j <- 0; 
				}			
				release list_detectedObjs; // empty list of detected objs
			}
			else
			{
				write "no objs perceived!";
			}
			
			//verify bestTarget
			if(bestTarget != nil)
			 {
			 	bTargetFound <- true;
			 }
			else
			{
				write "no target found! " + list_detectedObjs;
			}
			
			if (bVarComputeExecTime) {t16 <- t16 + machine_time - t;	}	
			
			return bTargetFound;
		}
		
		/**************************************************************************************************
		 * A7: fly randomly in the direction of the target with the amplitude of 30Â°left and 30Â° right
		 *************************************************************************************************/
		action flyToBestTarget
		{
		
			//check	target reached?
			if(bestTargetPoint distance_to self.location > TARGET_MAXDIST_TOLERENCE) //tolerance = 10 cm
			{
					heading <- self direction_to bestTargetPoint;
					point prv_loc <- location;
					//do update_limit;
					do wander speed: speed/step amplitude: 60 bounds: mySpatObj.geom_limit;
					do saveFlightDistance(prv_loc);
					iCmpTimeflyTo <- iCmpTimeflyTo + 1;					
			}
			else//reached it!
			{
				bBestTargetReached <- true;
				if(bCheckActStat)
				{
					do save_activity("flyTo",(iCmpTimeflyTo * step)); 
				}
				iCmpTimeflyTo <-0;
			}				
			fEnergyGainInAct <- 0.0;
			fEnergyLostInAct <- step * WANDER_ENERGLOST_RATE;					
			do updateCurrEnergy;
			do updateTargetWeight;	
		}
		
		/**********************************************************
		 * PRIVATE accessory METHODE
		 * calculate current energy from these sources:
		 * - stock of blood (fQntStckBlood)
		 * - stock of nectar (estimated from fTimeNeedAct to take nectar) 
		 * - time of resting 
		 * Energy unity is sec!
		 **********************************************************/	
		action updateCurrEnergy
		{
			fCurrEnergy <- (fCurrEnergy + fEnergyGainInAct) - fEnergyLostInAct;
			
			if(fCurrEnergy > MAX_ENERGY) 
			{
				fCurrEnergy <- MAX_ENERGY;
			}
			else if (fCurrEnergy <= 0)
			{
				do killMe("energy");
			}
			
		}
		
		/************************************************************
		 * calculate the effect of blood quatity to a target need
		 *************************************************************/
		float calStckBloodEff(string strTarget)
		{
			float fValue <- 0.0;
			if(fQntStckBlood > MAX_STOCK_BLOOD)
			{
				fQntStckBlood <- MAX_STOCK_BLOOD;
			}
			else if (fQntStckBlood < 0.00)
			{
				fQntStckBlood <- 0.0;
			}
			
			///
			switch(strTarget)
			{
				match "nectar"
				{
					float a <- 1.0;
					if(strState = "virgin")
					{
						float a <- 1.7;
					}
					fValue <- a*(0.5 - 0.3 * (fQntStckBlood/MAX_STOCK_BLOOD));
				}
				match "blood" 
				{
					fValue <- 1 - (fQntStckBlood/MAX_STOCK_BLOOD);
				}
				match "shade" 
				{
					fValue <- (fQntStckBlood/MAX_STOCK_BLOOD);
				}
				match "bs" 
				{
					if(strState = "gono")
					{
						fValue <-(fQntStckBlood/MAX_STOCK_BLOOD);
					} 
					else
					{
						fValue <- 0.0;
					}
				}
			}
			return fValue;
		}
		
		
		/**************************************************************
		 * calculate the effect of energy to a target need
		 *************************************************************/
		float calEnergyEff(string strTarget)
		{
			float fValue <- 0.0;
			switch(strTarget)
			{
				match "nectar"
				{
					fValue <- 0.8 - 0.8 * (fCurrEnergy/MAX_ENERGY);
				}
				match "blood" 
				{
					fValue <- 1 - (fCurrEnergy/MAX_ENERGY);
				}
				match "shade" 
				{
					fValue <- 0.483 - 0.3*(fCurrEnergy/MAX_ENERGY);
				}
				match "bs" 
				{
					if((fCurrEnergy/MAX_ENERGY)>= 0.5)
					{
						fValue <- (1.76*(fCurrEnergy/MAX_ENERGY)) - 0.818;
					} 
					else
					{
						fValue <- 0.0;
					}
				}

			}
			return fValue;
		}
		
		/*******************************************************
		 * PRIVATE ACTION: update target's desired weight 
		 ********************************************************/
		action updateTargetWeight
		{				
				//list_finalTargertNeedWg[mapTargetIndex["nectar"]] <- (BLOODEFF_WG_4_NECTAR * calStckBloodEff("nectar")) +((1-BLOODEFF_WG_4_NECTAR) + calEnergyEff("nectar"));
				list_finalTargertNeedWg[mapTargetIndex["nectar"]] <- calStckBloodEff("nectar")* calEnergyEff("nectar");
				list_finalTargertNeedWg[mapTargetIndex["blood"]] <- (BLOODEFF_WG_4_BLOOD*calStckBloodEff("blood")) + ((1-BLOODEFF_WG_4_BLOOD) * calEnergyEff("blood"));
				list_finalTargertNeedWg[mapTargetIndex["shade"]] <- (BLOODEFF_WG_4_SHADE*calStckBloodEff("shade")) + ((1-BLOODEFF_WG_4_SHADE) * calEnergyEff("shade"));
				//list_finalTargertNeedWg[mapTargetIndex["bs"]] <- (BLOODEFF_WG_4_BS * calStckBloodEff("bs")) + ((1-  BLOODEFF_WG_4_BS)*calEnergyEff("bs"));
				list_finalTargertNeedWg[mapTargetIndex["bs"]] <- calStckBloodEff("bs")*calEnergyEff("bs");
		}
		

		
		/**********************************************************************
		 * check if the physiodev for a gonotrophic female has already reached.
		 **********************************************************************/
		action checkPhysioDevDone
		{
		
			if(avgHourAirTempK_matrix[iSimRealDay, int(curTime)] >= MIN_TEMP_LAYEGG+ 273.15)//Â°K
			{ 	
				 if(iCmpGonoCycle > 1) //if the gono cycle is not the first time
				 {
				 	if(fCurrGonoPhysioDev >= LATE_GONO_THRES)
				 	{
				 		bGonoBegan <- false; //reinitialise the variable for the next gono cycle
				 		bPhysioDevDone <- true;
				 		fCurrGonoPhysioDev <- 0.0;
				 	}
				 }
				 else if (iCmpGonoCycle = 1)
				 {
				 	if(fCurrGonoPhysioDev >= FRS_GONO_THRES)
				 	{
				 		bGonoBegan <- false; //reinitialise the variable for the next gono cycle
				 		fCurrGonoPhysioDev <- 0.0;
				 		bPhysioDevDone <- true;
				 	}
				 }	 
			}
		}
		
		
		/**************************************
		 *  calculate the duration of each activity
		 **************************************/
		float getTimeNeed
		{
			switch strTarget
			{
				match "nectar" {return rnd((MAX_NEEDTIME_NECTAR - MIN_NEEDTIME_NECTAR)* 1000)/1000 + MIN_NEEDTIME_NECTAR;}
				match "blood" {return rnd((MAX_NEEDTIME_BLOOD - MIN_NEEDTIME_BLOOD)* 1000)/1000 + MIN_NEEDTIME_BLOOD;}
				match "bs" {return rnd((MAX_NEEDTIME_LAYEGG - MIN_NEEDTIME_LAYEGG)* 1000)/1000 + MIN_NEEDTIME_LAYEGG;}
			}
		}
		
		
		/************************************************
		* Give target attemp probability !! To be added
		*************************************************/
		bool applyAttempSurvivalRate
		{
			switch strTarget
			{
				match "blood" {return flip(BITE_ATTEMP_SURV);}
				default {return true;}
			}
		}
		
		
		/***************************************************************************************
		 * save duration of each activity for each mosquito
		 * One mosquito for one map_actTotDailyDur. THis variable will be initiate each day by the global statement
		 ***************************************************************************************************************/
		action save_activity(string activity, float duration){

			iCmpAct <- iCmpAct + 1;
			//save total duration per action in order to make a daily statistic
			map_actTotDailyDur[activity] <- map_actTotDailyDur[activity]  + duration;//result = [flyTo,180.0; blood,42.759; ...]; 
			dailyAct << [activity, duration, hour_of_a_day, mn_of_an_hour, location, name, iCmpAct];
			
			switch(activity)
			{
				match "blood" { iCmpDailyTb <-iCmpDailyTb +1; }
				match "nectar" {iCmpDailyTn <-iCmpDailyTn +1;}
				match "bs" {iCmpDailyBs <-iCmpDailyBs +1;}
				match "shade" {iCmpDailyShd <-iCmpDailyShd +1;}
				match "flyTo" {iCmpDailyFt <-iCmpDailyFt +1;}
				match "wander" 
				{
					iCmpDailyWd <-iCmpDailyWd +1;
				}
			}
		}
		
		
		/****************************************************************
		 * save a maximum flight distance of each Aedes if ask for
		 *************************************************************/
		action saveFlightDistance(point prv_loc)
		{
			float dist <- prv_loc distance_to location;
			float dispersDist <- bornLoc distance_to location;
			
			fDailyFlightDist  <- fDailyFlightDist + dist;
			if(dispersDist > fMaxLifeDisps)
			{
				fMaxLifeDisps <- dispersDist;
				fCoordX <- location.x;
				fCoordY <- location.y;
			}
			
			if(bTrackAedes)
			{
				add prv_loc to: movingPts_list;
				add location to: movingPts_list;			
			}
		}
		
		//**************************************************
		//    appearance of mosquito in gerenal
		//**************************************************
		aspect default
		{

			if(strState = "virgin")
			{
				draw circle(0.5) color: rgb("blue");
			}
			else if(strState = "ovipo")
			{
				draw circle(0.5) color: rgb("pink");
			}
			else if (strState = "gono")
			
			{
				draw circle(0.5) color: rgb("magenta");
			}

			//------- draw a moving trail -------
			if(bTrackAedes) 
			{
	 			draw trail color: rgb("blue");
	 			draw perception_area  color: rgb("cyan");
	 			draw text: string(strTarget) size: 20 color: rgb('black'); 
	 		}			
		}
	}
	
	
	/*******************************************************
 	* 				SPECIES SpatObj
 	******************************************************/
 	/////////////// parent agent for all kind of land use /////////////////
	species SpatObj
	{
		geometry geom_limit;	
		
		Dispersal myMaxDisp;		
		string OBJID;
		Class_Landuse spcClass;
		string strClass;
		
		int iSpatObjAirMinTempC; //use integer to decrease the sensibility of temperature comparation in decimal
		int iSpatObjAirMaxTempC;
		float fSpatObjAirAvgTempC;
		
		float fSpatObjWaterTempCIn;
		float fSpatObjWaterTempCOut;
		
		int iNbBsInNoWaterMax;
		int iNbBsInNoWaterMin;
		int iNbBsOutNoWaterMax;
		int iNbBsOutNoWaterMin;
		
		//-------------  for bBsVarMode only -------
		int iNbBsInNoWater;
		int iNbBsOutNoWater;
		//-------------------------------------------------
	
		int iNbBsInWithWater;
		int iNbBsOutWithWater;
			
		
		int iCapaHaveEgg;
		int iWaitStckEggIn <- 0;
		int iWaitStckEggOut <- 0;
		int iTotWaitStck update: iWaitStckEggIn + iWaitStckEggIn;
		list<float> fListStkEtoLIn <- []; //list's size 8 (maximum days to embryonate without death = 13Â°C)
		list<float> fListStkLtoPIn <- []; // list's size = 30? -at 13Â°C
		list<float> fListStkPtoAIn <- []; // list's size = 7 at 13Â°C
		
		list<float> fListStkEtoLOut <- []; //list's size 8 (maximum days to embryonate without death = 13Â°C)
		list<float> fListStkLtoPOut <- []; // list's size = 30? -at 13Â°C
		list<float> fListStkPtoAOut <- []; // list's size = 7 at 13Â°C
		
		int dayNeedEmbryo; //number of days needed to transit from eggs - larvae
		int dayNeedPupat; // from larvae - pupae
		int dayNeedEmerg; // from puape - adults

		float fQntHatchEggRateIn;
		float fQntHatchEggRateOut;
		float fEggSurvTempRate;
		float fLarvaSurvTempRate;
		float fPupaSurvTempRate;
		float fAdultSurvTempRate;
		int iNbNewFemaleA; //total
		int iNbNewFemaleAIn;
		int iNbNewFemaleAOut;

		int C2 <- 2; //must be calibrated such that the hatching increases more rapidely for highter value of wd, thus will create a quadratic-like function
		float WDMAX <- 70.0;	//for Delhi
      	
		float fSunExpoIn;
		float fSunExpoOut;
      	
		float fNbNewE;
		float fNbNewL;
		float fNbNewP;
		
		float fTotStckE;
		float fTotStckP;
		float fTotStckL;
		float totPopAquatic;
		
		list<float> list_attrctTarget<- [1.0,1.0,1.0,1.0,1.0];// ref mapTargetIndex["nectar"::0, "blood"::1, "shade"::2, "bs"::3, "not"::4];
  
 		list<SpatObj> my_neighbors; //returns all the agents located at a distance lower or equal to 10 to the agent applying the operator considering
		
		float fInitPop;
		float fCurrNbResidPop; // % of pop presence in a residential surface of the spatial object
		float fCurrNbCommPop; // % of pop presence in a commercial (pulbic) surface of the spatial object
		float fCurrNbTransPop; //% of population presence in a transportation surface of the spatial object
		
		int iLightProfile;  	
  		rgb color;
  		float fPorosity;
  		float fMaxDist;
  		
  		//--------mosquito activity check --------
  		int iHourBite; // increasing values by aedes
  		int iHourLayEgg;
  		int iTotBites;
  		int iTotLayEggs;
  		float fHourDensPop; //density of human per an hour
  		
  		reflex saveMosAct when: bCheckNbBiteOvipoPerHourInObj and bNewHour
  		{
  		  	float t <- bVarComputeExecTime ? machine_time : 0;  
  			
  			iTotBites <-iTotBites + iHourBite;
  			iTotLayEggs <- iTotLayEggs + iHourLayEgg;
  			
  			if(iHourBite != 0 or iHourLayEgg!= 0 )
  			{
	  			save[name,nb_days, hour_of_a_day, iHourBite, iHourLayEgg] type: csv to: "NbBiteOvipoInSite.csv";
	  			iHourBite <- 0;
	  		 	iHourLayEgg <- 0;
	  		}
	  	  	if (bVarComputeExecTime) {t7 <- t7 + machine_time - t;  }
  		}

		//-------------------------------------------------------
		// update the dispersal circle of a mosquito
		//----------------------------------------------------
		action updateDispersalCircle(Aedes an_aedes)
		{
			fMaxDist <- max([fMaxDist, an_aedes.fMaxLifeDisps]);
			myMaxDisp.fMaxRadius <- fMaxDist;
		}
  		
  		//-------------------------------
		// intanciate SpatObj agent
		//-----------------------------
		action initialisation
		{		 
			my_neighbors <- SpatObj at_distance 1.0;//topology(self) neighbours_of self; //returns all the agents located at a distance lower or equal to 10 to the agent applying the operator considering			
			add self to: my_neighbors;
			geom_limit <- union(my_neighbors);
			
			ask world{myself.color <-  adjustObjColor(myself.spcClass.color);}
			do calDailySpatObjAirTempC;
			do calDailySpatObjWaterTempC;
			
			//calculate the absolute number of breeding site from the density found in land use class of a spatial object, multiply by its surface
			if(bAvgDenseBSMode)
			{
				iNbBsInNoWater <- int(spcClass.fDenBsInNoWater * shape.area); //hasard sur la densité moyenne
				iNbBsOutNoWater <- int(spcClass.fDenBsOutNoWater * shape.area);
			}
			else
			{
				do calNbBsNoWater;			
			}
			
			fInitPop <- spcClass.fRecCapaRate * shape.area;
			fHourDensPop <- fInitPop/shape.area;			
			do updateBloodAttract;
			do updateNectarAttract; // static target
			do updateRestAttract;
			do updateBsAttract;
							
			dayNeedEmbryo <- int(rnd((MAX_DAYNEED_EMBRYONATION - MIN_DAYNEED_EMBRYONATION) * 1000)/1000 + MIN_DAYNEED_EMBRYONATION); //number of days needed to transit from eggs - larvae
			dayNeedPupat<- int(rnd((MAX_DAYNEED_PUPATION - MIN_DAYNEED_PUPATION) * 1000)/1000 + MIN_DAYNEED_PUPATION);  // from larvae - pupae
			dayNeedEmerg <- int(rnd((MAX_DAYNEED_EMERGENCE - MIN_DAYNEED_EMERGENCE) * 1000)/1000 + MIN_DAYNEED_EMERGENCE); 
			
			//fFloodHatchRatio <- rnd(EGG_FLOOD_HATCH_RATIO_MAX - EGG_FLOOD_HATCH_RATIO_MIN * 1000)/1000 + EGG_FLOOD_HATCH_RATIO_MIN; //FLOOD_HATCH_RATIO value uniform random (CimSim)
			
			// initialize the list of eggs to embryonate
			loop i from: 0 to: MAX_DAYNEED_EMBRYONATION-1
			{
				add 0.0 to: fListStkEtoLIn;
				add 0.0 to: fListStkEtoLOut;
			}
			
			// initialize the list of larvae to pupate
			loop i from: 0 to: MAX_DAYNEED_PUPATION-1
			{
				add 0.0 to: fListStkLtoPIn;
				add 0.0 to: fListStkLtoPOut;
			}	
					
			// initialize the list of larvae to pupate
			loop i from: 0 to: MAX_DAYNEED_EMERGENCE-1
			{
				add 0.0 to: fListStkPtoAIn;
				add 0.0 to: fListStkPtoAOut;
			}
			
			fTotStckE <- 0.0;
			fTotStckP <- 0.0;
			fTotStckL <- 0.0;
			totPopAquatic <- 0.0;
			
			if(!bOpenMole)
			{
				fPorosity <- rnd(spcClass.fPorosityMax * 1000)/1000 + spcClass.fPorosityMin;
			}
			else
			{
				fPorosity <- spcClass.fPorosityAvg;
			}		
		}
		
		
		//------------------------------------------------------------------------
		//  REFLEXES
		//  dynamic targets which need to update everyday (resting zone and bs)
		//------------------------------------------------------------------------
		reflex updateDynmTargetsAttrc when: bNewDay
		{
			float t <- bVarComputeExecTime ? machine_time : 0;  
			do updateRestAttract;
			do updateBsAttract;	
			if (bVarComputeExecTime) {t9 <- t9 + machine_time - t;}
		}
				
		/**********************************
		 * A NORMALISER LA VALEUR
		 **********************************/
		reflex updateBloodAttract when: bNewHour
		{
			float t <- bVarComputeExecTime ? machine_time : 0;  
			do updateBloodAttract; 		
			if (bVarComputeExecTime) {t10 <- t10 + machine_time - t;	}		 			 		
	 	}		
			
		/********************************************
		 * update percentage of Aedes actif by hour
		 * 
		 *******************************************/
		reflex checkAedesActifPeriod when: bNewHour
		{
			float t <- bVarComputeExecTime ? machine_time : 0; 
			float fActifRate <- 100.0;
			float fInRate<- 1.0;
			float fSocioRate <- 1.0;
			int iNextHour <- hour_of_a_day + 1;
			
			//to avoid error of running out of bound
			if(hour_of_a_day = 23)
			{
			 	iNextHour <- 0;	//avoid error of size 24!
			}
			
			//for the peak period
			fActifRate <-list_AeActifRateByHour[hour_of_a_day];
			
		    //day resting period (2h after sunrise and 2h before sunset) -> %Aedes actif can be actif at inside surface
			if((curTime > ((list_sunrise at (nb_days)) + 2)) and (curTime < ((list_sunset at (nb_days)) - 2)))
			{		
				// apply function regarding to the location inside or outside
				if(spcClass.fOutSpaceRate < 60)//inside
				{
					fInRate <-  (((-1.153846154) * (hour_of_a_day ^2) + 27.85714286*hour_of_a_day - 134.2857143)/100)+0.6; // use polynominial function
					fActifRate <- 0.8*(list_AeActifRateByHour[hour_of_a_day] + (((list_AeActifRateByHour[iNextHour]-list_AeActifRateByHour[hour_of_a_day])/60)* mn_of_an_hour)) + 0.2*fInRate;					
					list_attrctTarget[mapTargetIndex["shade"]]  <- 1.0;
				}
			}	
			//night resting (this method will work till midnight 
			else if ((curTime > (list_sunset at (nb_days))) and (curTime < 24)) // ((curTime < iNextDaySunriseHour))) //lay egg must do after sunset
			{
				if(spcClass.socialClass	= "low") //assume that lower classe doen't have electricity, the rate is slower than usually
				
				{
					fSocioRate <-0.1* (24-hour_of_a_day); 
					fActifRate <- list_AeActifRateByHour[hour_of_a_day]*fSocioRate;	
					list_attrctTarget[mapTargetIndex["shade"]]  <- 1.0; //change space to resting zone
				}
			}
			else if(curTime < iNextDaySunriseHour) //midnight  1am, 2am, ... till next sunrise
			{
				list_attrctTarget[mapTargetIndex["shade"]]  <- 1.0; //change space to resting zone	
			}

			//stochasetic
			fActifRate <- rnd(fActifRate*1000)/1000 + 0.05; 
			
			// wake up Aedes
			ask  Aedes overlapping self
			{
				bool bWakeup <- flip(fActifRate);
				
				 if(bWakeup and !bActive)
				 {
					 do wakeUp;	
				 }	
				 else if (bActive and !bWakeup)	//make rest for the active aedes
				 {
						list_finalTargertNeedWg[mapTargetIndex["shade"]] <- 1.0; //rising the priority to rest
				 }
			}	
			if (bVarComputeExecTime) {t11 <- t11 + machine_time - t;}		

		}
				
		
		/******************************************************************************
		 * REFLEX: Update stock of Aedes in aquatic stages and creat new female adulte
		 ******************************************************************************/
		reflex updateStocksAedes when: bUpdateStck and (bNewDay and ((iWaitStckEggIn != 0) or (iWaitStckEggOut != 0) or (totPopAquatic != 0.0)))
		{
			float t <- bVarComputeExecTime ? machine_time : 0;  
			fNbNewE <- 0.0;
			fNbNewL <- 0.0;
			fNbNewP <- 0.0;
			iNbNewFemaleA <- 0;
			
			// for memory 
			float fSaveNewE <- 0.0; 
			float fSaveNewL <- 0.0;
			float fSaveNewP <- 0.0;
						
			do calDailySpatObjAirTempC;
			do calDailySpatObjWaterTempC;			
			do calHatchEggQntRate; //give value to fQntEggHatchIn/Out
      		
			//------------------- stock of stage development --------------------
			bool bGoodTempDevIn  <- getStageDevDuration(fSpatObjWaterTempCIn, true);
			bool bGoodTempDevOut <- getStageDevDuration(fSpatObjWaterTempCOut, false);
			
			//===================
			// ------- Inside ------
			//===================			
			if(bGoodTempDevIn) //temperature dependent function check, if the temperature is so cold or so hot, no dev occurs
			{
				do checkTemperatureAquaticSurvivalRate(true); //calculate survival rate for eggs, larvae, pupae for inside water temp
				fSaveNewE <- updateStckEtoL(fListStkEtoLIn, iWaitStckEggIn, true);
				iWaitStckEggIn <- int(iWaitStckEggIn - fSaveNewE); //decrease the stock of waiting eggs from those who are submerged.
				fSaveNewL <- updateStckLtoP(fListStkLtoPIn,fListStkEtoLIn);
				fSaveNewP <- updateStckPtoA(fListStkPtoAIn, fListStkLtoPIn);	 			
			}
			
			iNbNewFemaleAIn <- int(fListStkPtoAIn[0]* EMERGC_SURV * FAMALE_PROB);		
			if(bCheckStck)//if save to files
			{
				string tmp <- "" + nb_days +","+ name + ","+ strClass+ ","+ fSaveNewE+ ","+ fSaveNewL+ ","+ fSaveNewP+ "," + iNbNewFemaleAIn+ "," +  fSpatObjWaterTempCIn + ", in" + "\n";
				txtStock <-txtStock + tmp; 
			}			
			
				
			//=========================
			//------- Outside -------
			//=========================
			if(bGoodTempDevOut) //temperature dependent function check, if the temperature is so cold or so hot, no dev occurs
			{		
				do checkTemperatureAquaticSurvivalRate(false); //assign survival rate for eggs, larvae, pupae for outside water temperature
				fSaveNewE <- updateStckEtoL(fListStkEtoLOut, iWaitStckEggOut, false);
				iWaitStckEggOut <- int(iWaitStckEggOut - fSaveNewE);
				fSaveNewL <-  updateStckLtoP(fListStkLtoPOut, fListStkEtoLOut);
				fSaveNewP <- updateStckPtoA(fListStkPtoAOut, fListStkLtoPOut);				
			}
			
			iNbNewFemaleAOut <- int(fListStkPtoAOut[0] * EMERGC_SURV * FAMALE_PROB);		
			if(bCheckStck)//if save to files
			{
				string tmp <- "" + nb_days +","+ name + ","+ strClass+ ","+ fSaveNewE+ ","+ fSaveNewL+ ","+ fSaveNewP+ "," + iNbNewFemaleAOut+ "," +  fSpatObjWaterTempCOut + ", out" + "\n";
				txtStock <-txtStock + tmp; 					
			}
							
			iNbNewFemaleA <- iNbNewFemaleAIn + iNbNewFemaleAOut;
			//int((fListStkPtoAIn[0] + fListStkPtoAOut[0]) * EMERGC_SURV * FAMALE_PROB);		
			fListStkPtoAIn[0] <- 0.0; //assume all pupae is emerged, the different of newFemal and the stock is males
			fListStkPtoAOut[0] <- 0.0; //assume all pupae is emerged, the different of newFemal and the stock is males
			fTotStckE <- sum(fListStkEtoLIn) + sum(fListStkEtoLOut);
			fTotStckL <- sum(fListStkLtoPIn) + sum(fListStkLtoPOut);
			fTotStckP <- sum(fListStkPtoAIn) + sum(fListStkPtoAOut);
			totPopAquatic <- fTotStckE + fTotStckL + fTotStckP;	


			// create a new agent AEDES - COMMENT this function if need the simulation run only for ONE GENERATION
			create Aedes number: iNbNewFemaleA
			{
				location <- any_location_in(myself); //any_location_in(initialBS);
				bornLoc <- location;
				iAge <- 0;	
				bornObj <- myself;
			}
			
			//if dispersal range is needed to be saved, check if the dispersal agent already exist, if not create a new one
			if(bDisprsShpNeed and myMaxDisp = nil)
			{
					create Dispersal 
					{
						strOrigSpatObj <- name; 
						myself.myMaxDisp <- self; // assign this dispersal agent to spatobj zho created the dispersal
					}
			}
						
			do killAdultAedesWithTempC;
			if (bVarComputeExecTime) {t12 <- t12 + machine_time - t;	}
		}
				
		/************************************
		 * update the attraction for nectar
		 ************************************/
		action updateNectarAttract
		{
			list_attrctTarget[(mapTargetIndex["nectar"])] <- spcClass.fNectarSourceRate;	//A rÃ©cupÃ©rer dans le shapefile via la valeur NDVI, pour l'instant rÃ©cupÃ©rÃ© via la classe
		}
		
		
		/******************************************************************************
		 * update blood attraction as a function of social class and daily temperature
		 *****************************************************************************/
		action updateBloodAttract
		{
		
			if(fInitPop > 0.000000) //if the object is not a vegetation
			{
				fCurrNbResidPop <- spcClass.fsurfPrivRate * fInitPop * calPopResidVarRate();
				fCurrNbCommPop <-	spcClass.fSurfPubRate * fInitPop * calPopCommVarRate();
				fCurrNbTransPop <- spcClass.fSurfTransRate * fInitPop * calPopTransVarRate();
				fHourDensPop  <- (fCurrNbResidPop + fCurrNbCommPop + fCurrNbTransPop)/shape.area;
				list_attrctTarget[mapTargetIndex["blood"]]  <-fHourDensPop; 
	 		}
	 		else
	 		{
	 			list_attrctTarget[mapTargetIndex["blood"]] <- 0.0; 
	 		}		
		}
		
		
		/*******************************************************************************************************************
		 * Human population variation depending on residential surface of a spatial object and time
		 ******************************************************************************************************************/
		float calPopResidVarRate
		{
			return list_percentPresPopResid[int(curTime)];
		}
		
		/********************************************************************************************************************
		 * Human population variation depending on public surface of a spatial object and time
		 ********************************************************************************************************************/
		float calPopCommVarRate
		{
			return list_percentPresPopComm[int(curTime)];
		}
		
		
		/**********************************************************************************************
		 * Human population variation depending on transport surface of a spatial object and time
		 *********************************************************************************************/
		float calPopTransVarRate
		{
			return list_percentPresPopTrans[int(curTime)];
		}
			
				
		/***************************************************************************************************
		 *  The BS attraction level depend on the proportion of water-filling BS and the total number of potential BS
		 ***************************************************************************************************/
		action updateBsAttract //  = cal_iNbBsInWithWater
		{
				float fRateVarBsIn <- 0.0; // depending to fGlobalTodayTempC
				float fRateVarBsOut <- 0.0; // depending to fGlobalTodayRainfall
				
				float fDepRain <- 0.0; // for outside bs
				float fDepSocio <- 0.0; //for inside bs
				float fHazard <- 0.05;
				
				int iInitBsIn <- iNbBsInNoWater; //for avg bs mode only
				int iInitBsOut <- iNbBsOutNoWater; // for avg bs mode only
				
				fDepSocio <- calEffSocioPracRate(fGlobalTodayAirTempC); //CALCULAT water storage INSIDE (depending on temperature)		
				fDepRain <- calRainEffRate();//cal rain effect for OUTSIDE
				
				
				if(!bAvgDenseBSMode)
				{
					iInitBsIn <- iNbBsInNoWaterMin;
					iInitBsOut <- iNbBsOutNoWaterMin;
					//BS WITH WATE INSIDE
					iNbBsInWithWater <- int(iInitBsIn * (1 + fDepSocio));// calculer Ã  partir de la valeur minimal (c'est pour Ã§a que je fasse + valeur min)	
					if(iNbBsInWithWater > iNbBsInNoWaterMax){iNbBsInWithWater <- iNbBsInNoWaterMax;}
					else if (iNbBsInWithWater < iNbBsInNoWaterMin){iNbBsInWithWater <- iNbBsInNoWaterMin;}
					//BS WITH WATER OUTSIDE
					iNbBsOutWithWater <- int(iInitBsOut * (1 + fDepRain));// quantitÃ© minimale + nb calculÃ© Ã  base de quant min
					if(iNbBsOutWithWater > iNbBsOutNoWaterMax){iNbBsOutWithWater <- iNbBsOutNoWaterMax;}
					else if (iNbBsOutWithWater < iNbBsOutNoWaterMin){iNbBsOutWithWater <- iNbBsOutNoWaterMin;}
					//cal atraction based on filling bs/all maximum bs
					int iMaxAllBs <- iNbBsOutNoWaterMax + iNbBsInNoWaterMax;
					list_attrctTarget[mapTargetIndex["bs"]] <- iMaxAllBs = 0 ? 0 : (iNbBsOutWithWater + iNbBsInWithWater)/iMaxAllBs;	
				}
				else//avg bs mode
				{
					iNbBsInWithWater <- int(iInitBsIn * fDepSocio);//BS WITH WATE INSIDE
					iNbBsOutWithWater <- int(iInitBsOut * fDepRain);//BS WITH WATER OUTSIDE
					
					int totInitBsNoWater <-(iInitBsIn + iInitBsOut); 
					if(totInitBsNoWater > 0) 
					{
						list_attrctTarget[mapTargetIndex["bs"]] <- (iNbBsInWithWater + iNbBsOutWithWater)/totInitBsNoWater;	
						if(list_attrctTarget[mapTargetIndex["bs"]]  > 1.0){list_attrctTarget[mapTargetIndex["bs"]]  <- 1.0;}
					}else
					{
						list_attrctTarget[mapTargetIndex["bs"]] <- 0.0;
					}
				}	
		}
		
		
		/******************************************************************************************************
		 * attraction des zone de repos est liÃ©e au niveau de la lumiÃ¨re et la tempÃ©rature (moins chaud) de l'espace
		 *******************************************************************************************************/
		action updateRestAttract
		{
			//give an atrraction value for resting zone (shady)
			list_attrctTarget[mapTargetIndex["shade"]] <- spcClass.fRestSourceRate; // seul le % de prÃ©sence de zone ombragÃ©e suffit			
	
		}		
		

		/******************************************************************************************************************
		 * value for egg hatching probability (rainfall dependent) [equation adapted from Yusoff,2012 (calculation of G(t)]
		 ******************************************************************************************************************/
		float calRainEffRate
		{
			float fValue <- EGGS_HATCH_NO_FLOODING_PROB ; //// hatch if no rain or over flood in a BS a cause of heavy rainfall with EGGS_HATCH_NO_FLOODING_PROB 
			//correction de fonction pour qu'il correspond mieux à Yusoff - corrigé le 31/03/2015
			if(fGlobalTodayRainfall > 0) and (fGlobalTodayRainfall <= WDMAX)
			{
				fValue <- (fGlobalTodayRainfall / WDMAX)^C2;//c2=2, wdmax = 70mm
				if(fValue > 1.0) {fValue <- 1.0;} // impossible in this case
			}
			
			return fValue;
		}			
		
		/************************************************************************************
		 * Daily temperature can make a variation of BS.
		 * We can suppose that hotter temperature it is, more coolers are used, so more BS.
		 * This variation should depend on the land use class
		 * The methode use global air temperature by default (if there's no parameter passing)
		 *
		 **************************************************************************************/
		float calEffSocioPracRate(float fAirTemp <- fGlobalTodayAirTempC)
		{
				float fValue <- 0.2; //20% de chance pour les autres classes
				if((spcClass.socialClass = "low" or spcClass.socialClass = "middle") and fAirTemp >= MIN_TEMP_USECOOLER)
				{
					fValue <- (fAirTemp/MAX_WORLDTEMPC)^2; // a voir l'Ã©quation
				}
				return fValue;
		}
				
		
		/***********************************************************************************
		 * calculate the number of potentil BS from the density of BS attached to the class
		 ***********************************************************************************/
		action calNbBsNoWater
		{
			//make a random min density between density min & max
			float fDenBsInNoWaterMin <-  rnd(((spcClass.fDenBsInNoWaterMax - spcClass.fDenBsInNoWaterMin)/2)*100)/100 + spcClass.fDenBsInNoWaterMin;
			float fDenBsOutNoWaterMin <-  rnd(((spcClass.fDenBsOutNoWaterMax- spcClass.fDenBsOutNoWaterMin)/2)*100)/100 + spcClass.fDenBsOutNoWaterMin;
			
			iNbBsInNoWaterMax <- int(spcClass.fDenBsInNoWaterMax* shape.area);//+0.5
			iNbBsInNoWaterMin <- int(fDenBsInNoWaterMin* shape.area);//+0.5
			iNbBsOutNoWaterMax <- int(spcClass.fDenBsOutNoWaterMax* shape.area); //+0.5
			iNbBsOutNoWaterMin <- int(fDenBsOutNoWaterMin* shape.area);//+0.5
		}

		/********************************************************************************************************************
		 * Eggs can hatch when the BS is filled with water. We must distinguish the filling method according to the location
		 * of the BS (in = precipitation, out = temperature, MoMa)
		 * fct(fGlobalRainfall, fGlobalTemperarue, un paramÃ¨tre de vidange/remplissage journalier (random 100 <= 5)
		 **********************************************************************************************************************************/
		action calHatchEggQntRate
		{
			fQntHatchEggRateIn <-  max([EGGS_HATCH_NO_FLOODING_PROB,NOMINAL_DAILYSURV_EGGS_RATE * calEffSocioPracRate(0.0)]);//(0.5* EGGS_HATCH_NO_FLOODING_PROB) + 0.5 *(NOMINAL_DAILYSURV_EGGS_RATE * calEffSocioPracRate()) ;
			fQntHatchEggRateOut <- max([EGGS_HATCH_NO_FLOODING_PROB,NOMINAL_DAILYSURV_EGGS_RATE * calRainEffRate()]);//		 * value for egg hatching probability (rainfall dependent) [equation adapted from Yusoff,2012 (calculation of G(t)]	
		}
		
		/******************************************
		 * This fonction is called by AEDES agent
		 ******************************************/
		 float getCapaHaveEgg
		 {
		 	iCapaHaveEgg <- int(((iNbBsInWithWater + iNbBsOutWithWater)* AVG_EGG_BSSURFACE) - (iWaitStckEggIn + iWaitStckEggOut + sum(fListStkEtoLIn) + sum(fListStkLtoPIn) + sum(fListStkPtoAIn)+ sum(fListStkEtoLOut) + sum(fListStkLtoPOut) + sum(fListStkPtoAOut)));
		 	return iCapaHaveEgg;
		 }
		
		/********************************************************************************************
		 * Update stock of submerged eggs that can transform to lava stage at the end of their development 
		 *********************************************************************************************/
		float updateStckEtoL(list<float> fListStkEtoL, int iWaitStckEgg, bool bIn)
		{
				float fEggsToHatch <- 0.0;
				
				loop i from: 0  to: (length(fListStkEtoL) - 2)
				{
					fNbNewE<- NOMINAL_DAILYSURV_EGGS_RATE* fEggSurvTempRate* fListStkEtoL[i + 1];
					fListStkEtoL[i] <- (fListStkEtoL[i] + fNbNewE)  with_precision 1;
					fListStkEtoL[i+1] <- 0; //empty stock of the previous day
				}
				
				///add new stock
				if(iWaitStckEgg != 0)
				{
						if(bIn){
							fEggsToHatch <- fQntHatchEggRateIn * iWaitStckEgg;		
						} //if egg was laid inside, there's no effect of rain so fQntHatchEggRateOut = 1
						else
						{
							fEggsToHatch <- fQntHatchEggRateOut * iWaitStckEgg;			
						}
						fListStkEtoL[dayNeedEmbryo] <- (fListStkEtoL[dayNeedEmbryo] + fEggsToHatch);//  with_precision 0;
				}
				return fEggsToHatch;
		}
				
		/*******************
		 * updateStckLtoP 
		 ******************/
		 float updateStckLtoP(list<float> fListStkLtoP, list<float> fListStkEtoL)
		 {
		 		float addL <- 0.0;
		 		
		 		loop i from: 0  to: (length(fListStkLtoP) - 2)
				{
					fNbNewL<- NOMINAL_DAILYSURV_LARVAE_RATE* fLarvaSurvTempRate * fListStkLtoP[i + 1];
					fListStkLtoP[i] <- (fListStkLtoP[i] + fNbNewL) with_precision 1;
					fListStkLtoP[i+1] <- 0; //vider stock of the previous day (include mortalité)
				}
				
				//add the new larvae to the lastest postion on the list (depend on the day to pupate)	
				if(fListStkEtoL[0] != 0)
				{
					addL <- (EMBRYONATION_SURV* fListStkEtoL[0]);
					fListStkLtoP[dayNeedPupat] <- (fListStkLtoP[dayNeedPupat] + addL);//  with_precision 0;
					fListStkEtoL[0] <- 0;					
				}
				
				return addL;
		 }
	 	 
		 /***********************************
		  * update stock pupae to adulte 
		  ********************************/
		 float updateStckPtoA(list<float> fListStkPtoA, list<float> fListStkLtoP)
		 {
		 		float addP <- 0.0;
		 		loop i from: 0  to: (length(fListStkPtoA) - 2)
				{
					fNbNewP<- NOMINAL_DAILYSURV_PUPAE_RATE * fPupaSurvTempRate * fListStkPtoA[i + 1];
					fListStkPtoA[i] <- (fListStkPtoA[i] + fNbNewP)  with_precision 1;
					fListStkPtoA[i+1] <- 0; //vider stock of the previous day
				}
				
				if(fListStkLtoP[0] != 0)
				{
					addP <- (PUPATION_SURV* fListStkLtoP[0]);
					fListStkPtoA[dayNeedEmerg] <- (fListStkPtoA[dayNeedEmerg] + addP); //  with_precision 0;
					fListStkLtoP[0] <- 0;					
				}
				return addP;
		}	
		
		/******************************************************************************************************
		 * This action is called by the reflex update stock to kill some Aedes mosquitoes upon the temperature
		 *****************************************************************************************************/
		action killAdultAedesWithTempC
		{

	     	//---- ADULT : temperature dependency --------------
	     	if(fSpatObjAirAvgTempC > ADULT_DAILYSURV_LOW_TEMP_LIMIT and fSpatObjAirAvgTempC < ADULT_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT) // increase function
	     	{
	 			fAdultSurvTempRate <- calTempSurvRate (fSpatObjAirAvgTempC, ADULT_DAILYSURV_LOW_TEMP_LIMIT, ADULT_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT, ADULT_DAILYSURV_LOW_TEMP_RATE, 1.0);   		
	     	}
	     	else if (fSpatObjAirAvgTempC > ADULT_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT and fSpatObjAirAvgTempC < ADULT_DAILYSURV_HIGH_TEMP_LIMIT) //decrease fucntion
	     	{
	     		fAdultSurvTempRate <- calTempSurvRate (fSpatObjAirAvgTempC, ADULT_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT, ADULT_DAILYSURV_HIGH_TEMP_LIMIT, 1.0, ADULT_DAILYSURV_HIGH_TEMP_RATE);   		
	     	}
	     	else if (fSpatObjAirAvgTempC >= ADULT_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT and fSpatObjAirAvgTempC <= ADULT_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT)
	     	{
	     		fAdultSurvTempRate <- 1.0;
	     	}
	     	else if (fSpatObjAirAvgTempC <= ADULT_DAILYSURV_LOW_TEMP_LIMIT) 
	     	{
	     		fAdultSurvTempRate <- ADULT_DAILYSURV_LOW_TEMP_RATE;
	     	} 
	     	else if (fSpatObjAirAvgTempC >= ADULT_DAILYSURV_HIGH_TEMP_LIMIT)
	     	{
	     		fAdultSurvTempRate <- ADULT_DAILYSURV_HIGH_TEMP_RATE;
	     	}	     	
	     	
	     	int cmp <- 0;
	     	//kill some Aedes who are inside the SpatObj
	     	
	     	ask  Aedes overlapping self
	     	{
	     		if(!flip(myself.fAdultSurvTempRate))
	     		{
					do killMe("temperature");
					cmp <- cmp + 1; 
	     		}
	     	}
		}
				
		
				
		/*******************************************************************************************************************
		survival increases linearly with two extreme limits DAILY_MIN_TEMP from 0.05 at LOW_TEMP_LIMIT to 1.0 at HIGH_TEMP_LIMIT 
		(cf: Skeeter Buster (Text S2 : 
		Details of CIMSiM elements used in Skeeter Buster, together with modifications adopted) for the function
		C'est moi qui fait l'Ã©quation 
		*************************************************************************************************************************/
	    float calTempSurvRate (float temp, float temp1, float temp2, float rate1, float rate2)
	    {
	    	float a <- (rate2-rate1)/(temp2-temp1);
	    	float b <- rate1 - a*temp1;
			//write "temp surv " +  (a * temp) + b;
	    	return (a * temp) + b;
	    }
	    
		/********************************************************************************
		 * Survival rate of aquatic stages depends on SpatObj's water temerature
		 * @param: fSpatObjWaterTempC - temperature of water inside and outside
		 *******************************************************************************/
		action checkTemperatureAquaticSurvivalRate(bool bInside)
		{
	     	
	     	float fSpatObjWaterTempC <- 0.0 ;
			
			//assign a value
	     	if(bInside)
	     	{
	     		fSpatObjWaterTempC <- fSpatObjWaterTempCIn;
	     	}
	     	else
	     	{
	     		fSpatObjWaterTempC<- fSpatObjWaterTempCOut;
	     	}
	     	
	     	
			//----------- egg-----------------
	     	if(fSpatObjWaterTempC > EGG_DAILYSURV_LOW_TEMP_LIMIT and fSpatObjWaterTempC < EGG_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT) // increase function
	     	{
	 			fEggSurvTempRate <- calTempSurvRate (fSpatObjWaterTempC, EGG_DAILYSURV_LOW_TEMP_LIMIT, EGG_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT, EGG_DAILYSURV_LOW_TEMP_RATE, 1.0);   		
	     	}
	     	else if (fSpatObjWaterTempC > EGG_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT and fSpatObjWaterTempC < EGG_DAILY_SURVIVAL_HIGH_TEMP_LIMIT) //decrease fucntion
	     	{
	     		fEggSurvTempRate <- calTempSurvRate (fSpatObjWaterTempC, EGG_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT, EGG_DAILY_SURVIVAL_HIGH_TEMP_LIMIT, 1.0, EGG_DAILYSURV_HIGH_TEMP_RATE);   		
	     	}
	     	else if (fSpatObjWaterTempC >= EGG_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT and fSpatObjWaterTempC <= EGG_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT)
	     	{
	     		fEggSurvTempRate <- 1.0;
	     	}
	     	else if (fSpatObjWaterTempC <= EGG_DAILYSURV_LOW_TEMP_LIMIT) 
	     	{
	     		fEggSurvTempRate <- EGG_DAILYSURV_LOW_TEMP_RATE;
	     	} 
	     	else if (fSpatObjWaterTempC >= EGG_DAILY_SURVIVAL_HIGH_TEMP_LIMIT)
	     	{
	     		fEggSurvTempRate <- EGG_DAILYSURV_HIGH_TEMP_RATE;
	     	}
	     	
//	     	//---------- larvae ------------------
	     	if(fSpatObjWaterTempC > LARVAE_DAILYSURV_LOW_TEMP_LIMIT and fSpatObjWaterTempC < LARVAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT) // increase function
	     	{
	 			fLarvaSurvTempRate <- calTempSurvRate (fSpatObjWaterTempC, LARVAE_DAILYSURV_LOW_TEMP_LIMIT, LARVAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT, LARVAE_DAILYSURV_LOW_TEMP_RATE, 1.0);   		
	     	}
	     	else if (fSpatObjWaterTempC > LARVAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT and fSpatObjWaterTempC < LARVAE_DAILYSURV_HIGH_TEMP_LIMIT) //decrease fucntion
	     	{
	     		fLarvaSurvTempRate <- calTempSurvRate (fSpatObjWaterTempC, LARVAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT, LARVAE_DAILYSURV_HIGH_TEMP_LIMIT, 1.0, LARVAE_DAILYSURV_HIGH_TEMP_RATE);   		
	     	}
	     	else if (fSpatObjWaterTempC >= LARVAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT and fSpatObjWaterTempC <= LARVAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT)
	     	{
	     		fLarvaSurvTempRate <- 1.0;
	     	}
	     	else if (fSpatObjWaterTempC <= LARVAE_DAILYSURV_LOW_TEMP_LIMIT) 
	     	{
	     		fLarvaSurvTempRate <- LARVAE_DAILYSURV_LOW_TEMP_RATE;
	     	} 
	     	else if (fSpatObjWaterTempC >= LARVAE_DAILYSURV_HIGH_TEMP_LIMIT)
	     	{
	     		fLarvaSurvTempRate <- LARVAE_DAILYSURV_HIGH_TEMP_RATE;
	     	}
			
//	     	//---------- pupae ----------------------
	     	if(fSpatObjWaterTempC > PUPAE_DAILYSURV_LOW_TEMP_LIMIT and fSpatObjWaterTempC < PUPAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT) // increase function
	     	{
	 			fPupaSurvTempRate <- calTempSurvRate (fSpatObjWaterTempC, PUPAE_DAILYSURV_LOW_TEMP_LIMIT, PUPAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT, PUPAE_DAILYSURV_LOW_TEMP_RATE, 1.0);   		
	     	}
	     	else if (fSpatObjWaterTempC > PUPAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT and fSpatObjWaterTempC < PUPAE_DAILYSURV_HIGH_TEMP_LIMIT) //decrease fucntion
	     	{
	     		fPupaSurvTempRate <- calTempSurvRate (fSpatObjWaterTempC, PUPAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT, PUPAE_DAILYSURV_HIGH_TEMP_LIMIT, 1.0, PUPAE_DAILYSURV_HIGH_TEMP_RATE);   		
	     	}
	     	else if (fSpatObjWaterTempC >= PUPAE_DAILYSURV_NORMAL_TEMP_LOWER_LIMIT and fSpatObjWaterTempC <= PUPAE_DAILYSURV_NORMAL_TEMP_UPPER_LIMIT)
	     	{
	     		fPupaSurvTempRate <- 1.0;
	     	}
	     	else if (fSpatObjWaterTempC <= PUPAE_DAILYSURV_LOW_TEMP_LIMIT) 
	     	{
	     		fPupaSurvTempRate <- PUPAE_DAILYSURV_LOW_TEMP_RATE;
	     	} 
	     	else if (fSpatObjWaterTempC >= PUPAE_DAILYSURV_HIGH_TEMP_LIMIT)
	     	{
	     		fPupaSurvTempRate <- PUPAE_DAILYSURV_HIGH_TEMP_RATE;
	     	}

		}
		
		
		/*********************************************************************************
		 * check the duration for each stage development by checking the temperature of the day
		 * It must be a water temperature
		 *********************************************************************************/
		bool getStageDevDuration(float fSpatObjWaterTempC, bool bInside)
		{
			bool bDie <- true;

			if((fSpatObjWaterTempC > 0.0) and (int(fSpatObjWaterTempC) < length (temperature_nbdays_dev_matrix column_at 1)-1))
			{
				dayNeedEmbryo <- int((temperature_nbdays_dev_matrix row_at fSpatObjWaterTempC) at 2);
				dayNeedPupat<-  int((temperature_nbdays_dev_matrix row_at fSpatObjWaterTempC) at 3);
				dayNeedEmerg <-  int((temperature_nbdays_dev_matrix row_at fSpatObjWaterTempC) at 4);
				bDie <- true;
	
			}	
			else //if the temperature is too cold < 0Â°c or too hot > 50Â°C, kill all stock
			{
				if(bInside)
				{
					release fListStkEtoLIn; //list's size 8 (maximum days to embryonate without death = 13Â°C)
					release fListStkLtoPIn; // list's size = 30? -at 13Â°C
					release fListStkPtoAIn;  
				}
				else //outside
				{
					release fListStkEtoLOut; //list's size 8 (maximum days to embryonate without death = 13Â°C)
					release fListStkLtoPOut; // list's size = 30? -at 13Â°C
					release fListStkPtoAOut;  
				}
				bDie <- false;
			}
			return bDie;
		}
		
		
		/*********************************************************************************************************
		 * Sun Expotion level is fluctuated depending on hour of the day ( polynominale fonction).
		 *We assume inside the building, there is no sun-exposure, so level = 0 
		 ********************************************************************************************************/	
		reflex dynamicSunExpo when: bNewHour
		{
			//a faire  varier selon l'heure
			if(hour_of_a_day >= list_sunrise[nb_days] and hour_of_a_day <= list_sunset[nb_days])
			{
				fSunExpoOut <- 0.000181387*(hour_of_a_day^4) - 0.00896242*(hour_of_a_day^3) + 0.141628*(hour_of_a_day^2) - 0.768125*hour_of_a_day + 1.40564;
			}
			else
			{
				fSunExpoOut <- 1.0;
			}
			
			fSunExpoIn <- 0.0;
		}	
		
		/*********************************************************************************************************
		* using water temperature equation to calculate max temperature and min temperature of water in container
	 	* It is sun exposure - temperature dependancy equation.
	 	* Equation ref: Text S2: details of CIMSiM used in Skeeter Buster
	 	*
	 	**********************************************************************************************************/		
		action calDailySpatObjWaterTempC
		{
			
			//cal water temp for inside using air temperature of SpatObj
			float fDailyWaterMaxTemp <- 15.03 + (0.27 * iSpatObjAirMinTempC) + (0.01 * ((iSpatObjAirMaxTempC)^2)) + (7.69*fSunExpoIn);
			float fDailyWaterMinTemp <- 5.02 + (0.81* iSpatObjAirMinTempC) + (0.001* ((iSpatObjAirMaxTempC)^2)) - (1.36*fSunExpoIn);
			fSpatObjWaterTempCIn <- (fDailyWaterMaxTemp + fDailyWaterMinTemp)/2;
			
			//water temperature for ouside breeding site, using overal air temp!!
			fDailyWaterMaxTemp <- 15.03 + (0.27 * list_daily_min_temp[nb_days]) + (0.01 * ((list_daily_max_temp[nb_days])^2)) + (7.69*fSunExpoOut);
			fDailyWaterMinTemp <- 5.02 + (0.81* list_daily_min_temp[nb_days]) + (0.001* ((list_daily_max_temp[nb_days])^2)) - (1.36*fSunExpoOut);
		 	fSpatObjWaterTempCOut <- (fDailyWaterMaxTemp + fDailyWaterMinTemp)/2; 
		}
		
		
		
		/************************************************************************
		 * calculate SpatObj's temperature by using a weight of SpatObj'type
		 ************************************************************************/
		action calDailySpatObjAirTempC
		{
			iSpatObjAirMinTempC <- int(list_daily_min_temp[nb_days]  * spcClass.fTempVarRate);
			iSpatObjAirMaxTempC <- int(list_daily_max_temp[nb_days]  * spcClass.fTempVarRate);
			fSpatObjAirAvgTempC <- (iSpatObjAirMinTempC + iSpatObjAirMaxTempC)/2;
		}	 	
		
		
		//------------------------------------------------
		// ASPECT
		//-------------------------------------------------
		
		aspect bloodAtt
		{
			rgb tmpColor <- rgb(255-(fHourDensPop*1000),0,0);
			draw shape empty: false color:tmpColor; //spcClass.color; // depth: fHeight;

		}
		
		aspect nectarAtt
		{
			rgb tmpColor <- rgb(0,255-(list_attrctTarget[(mapTargetIndex["nectar"])]*100),10);
			draw shape empty: false color:tmpColor; //spcClass.color; // depth: fHeight;
		}
		
		aspect bsAtt
		{
			rgb tmpColor <- rgb(0,0,255-(((iNbBsOutWithWater + iNbBsInWithWater)/shape.area)*100));
			draw shape empty: false color:tmpColor; //spcClass.color; // depth: fHeight;
			//draw text: string(((iNbBsOutWithWater + iNbBsInWithWater)/shape.area) with_precision 1 ) size: 5 color: rgb('white'); 
		}
		
		aspect stock
		{
			bool bEmpty <- true;
			rgb tmpColor <- color;
			
			if(totPopAquatic != 0.0)
			{
				tmpColor <- rgb(255-totPopAquatic,0,0);
				bEmpty		<- false;		
			}
			
			draw shape border: rgb("white") empty: bEmpty color:tmpColor; //spcClass.color; // depth: fHeight;
			draw text: string(int(totPopAquatic))  size: 3 color: rgb('white'); 
		}
		
		/**********************************************************************************************
		 * show nb of total bite and hourly bite per site only if bCheckNbBiteOvipoPerHourInObjet = true
		 ***********************************************************************************************/
		aspect showBite
		{
			if(bCheckNbBiteOvipoPerHourInObj)
			{
				draw shape empty: false color:color; //spcClass.color; // depth: fHeight;
				if(iTotBites != 0)
				{
					draw circle(iTotBites/100) at: location empty:false color:rgb("red");
					//draw text: string(int(iTotBites))  size: 3 color: rgb('white'); 
				}
				
				if(iHourBite != 0)
				{
					draw circle(iHourBite) at: location empty:false color:rgb("yellow") ;
					draw string(int(iHourBite)) at: location size: 3 color: rgb('black'); 
				}
			}
		}	
		
		/**
		 *  show dispersal
		 */
		aspect showDispersal
		{
			draw shape empty: false color:color;
			if(bDisprsShpNeed)
			{
				draw myMaxDisp.shape empty: false color: myMaxDisp.color;
			}
			else if(fMaxDist != 0.0)
			{
				draw circle(fMaxDist) at: location empty:false color:rgb("blue") ;
			}
			
		}		
	}
	
	
	/*******************************************************
 	* 	SPECIE Class_Landuse
 	* save different types of land use 
 	* Types must be loaded at the beginning of the model 
 	********************************************************/
	species Class_Landuse
	{
		string strClassName;
		string socialClass;
		float fRecCapaRate;
		float fsurfPrivRate;
		float fSurfPubRate;
		float fSurfTransRate;
		float fNectarSourceRate;
		float fRestSourceRate;
		float fOutSpaceRate;	
		float fTempVarRate;	
		float fDenBsOutNoWaterMin;
		float fDenBsOutNoWaterMax;
		float fDenBsInNoWaterMin;
		float fDenBsInNoWaterMax;
		float fPorosityMin;
		float fPorosityMax;
		float fPorosityAvg; // only for testing in OpenMole
		rgb color;
		
		////////////// Variables used only for bsAttVarMode ///////////
		float fDenBsInNoWater;
		float fDenBsOutNoWater;

		
	}
	
	/***********************
	 * STATISTIC
	 *********************/
	 species Stat
	 {
		list<geometry> activities_geom;
		list<rgb> activities_rgb;
		list<list> activities <- [];
		//string type_act <- nil;
		list<list> list_act;
		float x_max ;
		float y_max ;
		float emplacement_cur <- 0.0;
		float size <- 0.0;
		Aedes myAedes;
		
		
		//list<list> day_actGeom;	
		map<point, string> legends;
					
		init
		{
			x_max <- world.location.x * 2;
			y_max <- world.location.y * 2;
			myAedes <- one_of(Aedes);	
		}
		
		/**************************************************************************************
		* Save mosquitoes' activity function in order to produce a shapefile of activity chaine
		/**************************************************************************************/
		reflex create_act_geom when: cycle > 1 and bNewDay
		{
			activities_rgb <- [];
			activities_geom <- [];
			legends <- map([]);
			x_max <- world.location.x * 2;
		    y_max <- world.location.y * 2;
		    emplacement_cur <- 0.0;
		    
		    
			//add the last action before update dailyAct 
			ask myAedes
			{			
				if(fTimeStartResting != 0.0 and strTarget = "shade")
				{
						
						iCmpAct <- iCmpAct + 1;
						dailyAct << [strTarget, time - fTimeStartResting, hour_of_a_day, mn_of_an_hour, location, name, iCmpAct];
						fTimeStartResting <- time;
						if(debug){write "new day dailyAct "+dailyAct;}
				}
				myself.list_act <- dailyAct;
			}
			
			float tot <- sum(list_act collect float(each[1]));
			if(debug){write "tot dur = " + tot;}
			
			float totSize <- 0.0;
			float dur <- 0.0;
			bool haut <- true;
			loop act over: list_act
			{
				string strAct <- string(act[0]);
				activities_rgb << activities_color[strAct];
				int hour <- int(act[2]);
				int mn <- int(act[3]);
				point locAct <- point(act[4]);
				string aedesName <- string(act[5]);
				int iActSeq <- int(act[6]);
				
				size <- float(act[1]) / tot * x_max;
				if (size <= 0.1) 
				{
					size <- 0.1;
				}
				dur <- dur + float(act[1]);
				if (float(act[1]) > 300) {
					legends[{emplacement_cur, haut ? (world.location.y - y_max / 10.0) : (world.location.y + y_max / 10.0)}] <- string((dur / 3600) with_precision 1);
					haut <- not haut;
				}
				if(debug){write "size " + size;}
			
				create activity_stat with: [aedesName::aedesName, type::strAct, duration::float(act[1]), beginning::dur - float(act[1]) , color::rgb(activities_color[strAct]), sim_day::nb_days,hour::hour, mn:: mn, shape::rectangle(size,y_max / 10.0 ) translated_to {emplacement_cur + (size / 2.0) , world.location.y}];		
				create activity_loc  with: [actType::strAct, aedesName::aedesName, sim_day::nb_days,hour::hour, mn:: mn, shape::locAct, iActSeq::iActSeq ];
				save activity_loc type:"shp" to: "activityLoc1.shp" with:[aedesName::"AEDES", actType::"TYPE", sim_day::"DAY_SIM", hour::"HOUR", mn::"MN"];

				//activities_geom << rectangle(size,y_max / 10.0 ) translated_to {emplacement_cur + (size / 2.0) , world.location.y};				
				emplacement_cur <- emplacement_cur + size;
				totSize <- totSize + size;
			}

			ask myAedes
			{
					dailyAct <- []; //init the list for the new day
			}
		}
	}

	/*************************************
	 * 
	 ************************************/
	species activity_stat 
	{
	 	string type;
	 	float duration;
	 	float beginning;
	 	rgb color;
	 	int sim_day;
	 	int hour;
	 	int mn;
	 	string aedesName;
	 	
	 	aspect default 
	 	{
	 		draw shape color: color;
	 	}	
	 }
	

	species dailyTrail 
	{
		string aedesName;
		geometry trail;
	}
		
	/****************************************
	 * save activiy location of a mosquito
	 ***************************************/ 
	species activity_loc 
	{
	 	string actType;
	 	string aedesName;
	 	rgb color;
	 	int sim_day;
	 	int hour;
	 	int mn;
	 	int iActSeq; //the sequence of the activity happened in the mosquito's life
	 	
	 	aspect default 
	 	{
	 		draw shape color: color;
	 	}	
	 }
	 
	 
	 /*****************************************
	  * save maximum dispersal for each Aedes
	  ******************************************/
	 species dispersal_loc
	 {
	 	string aedesName;
	 	int sim_day;
	 	float coorX;
	 	float coorY;
	 }
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 				BATCH EXPERIMENTS
///////////////////////////////////////////////////////////////////////////////////////////////
experiment 'Run10sims' type: batch repeat: 10 keep_seed: true until: length(Aedes) = 0 
{
	int cmp <- 1;
	action _step_ 
	{
		save txtSaveDeath type: "csv" to: deathFileName+"-"+ cmp;	
		save txtStock type: "csv" to: fileDailyStock+"-"+cmp;
		save txtEggLaid type: "csv" to: fileEggLaid+"-"+cmp;
		cmp <- cmp + 1;
	}
}

//////////////////////////////////////////////////////////
// 		GRAPHICAL EXPERIMENTS
/////////////////////////////////////////////////////////
experiment fuzzyAedes45 type: gui 
{
	
	output 
	{
		
		//-----------------------------------------------------------
		// display the biting frequence in every spatial object
		//------------------------------------------------------------
		display Bites2D
		{
			species SpatObj aspect: showBite; // transparency: 0.5;
			species Dispersal aspect: default transparency: 0.5;
			species Aedes aspect: default transparency: 0.5; //female_fly;	
		}	

		//-------------------------------------
		// show different dispersal range
		//------------------------------------
		display showDispersal type: opengl
		{
			species SpatObj aspect: showDispersal refresh: false; // transparency: 0.5;
			species Aedes aspect: default transparency: 0.5; //female_fly;
		}	

		//----------------------------------------------------
		// Display of mosquitoes'activity dynamics
		//----------------------------------------------------
		display chart 
		{
			chart "time consuming for a blood meal" type: series size: {0.5,0.5} position: {0.0, 0.0} 
			{ 
				data "min" value: stat_takeBlood[0] color:rgb(250,128,114);
				data "max" value: stat_takeBlood[1] color:rgb(233,150,122);
				data "mean" value: stat_takeBlood[2] color:rgb("red");
				data "1er quartile" value: stat_takeBlood[3] color:rgb(205,92,92);
				data "median" value: stat_takeBlood[4] color:rgb(139,0,0);
				data "3er quartile" value: stat_takeBlood[5] color:rgb(178,34,34);
			}
			
			
			chart "time consuming for resting" type: series size: {0.5,0.5} position: {0.5, 0.0} 
			{ 
				data "min" value: stat_rest[0] color:rgb(250,128,114);
				data "max" value: stat_rest[1] color:rgb(233,150,122);
				data "mean" value: stat_rest[2] color:rgb("red");
				data "1er quartile" value: stat_rest[3] color:rgb(205,92,92);
				data "median" value: stat_rest[4] color:rgb(139,0,0);
				data "3er quartile" value: stat_rest[5] color:rgb(178,34,34);
			}
						
			chart "Average frequence of an activity/day/Aedes" type: series size: {0.5,0.5} position: {0.0, 0.5} 
			{
   				data "take blood" value: iNbAvgDailyTb color: rgb("red");
				data "take nectar" value: iNbAvgDailyTn color: rgb("green");
				data "lay egg" value: iNbAvgDailyBs color: rgb("blue");
				data "rest" value: iNbAvgDailyShd color: rgb("black");
				data "fly to" value: iNbAvgDailyFt color: rgb("cyan");
				data "wander" value: iNbAvgDailyWd color: rgb("yellow");
			}
			
			chart "average daily flight distance" type: series size: {0.5,0.5} position: {0.5, 0.5} 
			{
				data "distance" value: fAvgDailyFlyDist color: rgb("magenta");
			}
		}
		
		//------------------------------------
		// display stock dynamics
		//----------------------------------
		display Stock refresh_every: 1440
		{
			chart "stock" type: series size: {1,1} position: {0.0, 0.0} 
			{		
				data "eggs" value: totEggs color: rgb('blue');
				data "larvae" value: totLarvae color: rgb('yellow');
				data "pupae" value: totPupae color: rgb('green');
				data "new female" value: totNewAdults color: rgb('red');	
				data "all adults" value: length(Aedes) color: rgb('black');
 			}
		}
	}
}


	