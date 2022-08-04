#include "pnet_model.h"

void pnet_model::Phenology(veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int GrowthPhase)
{
	//
	// Phenlology calculations for the PnET ecosystem model.
	//

	double OldFolMass, GDDFolEff, delGDDFolEff, FolMassNew, spring_mod_select, fall_mod_select,
		GDD, CDD, t0, t0_chill, T_base, T_opt, T_min, T_max, Tc, C_ini, Rc, k, Rf,
		F_crit1, F_crit2, C_req, fall_t0, fall_T_base, fall_F_crit1, fall_F_crit2, f_fixed1, f_fixed2;



	if (GrowthPhase == 1)
	{
		switch ( clim->timestep )
		{
			case 0:   // monthly
				share->dayspan = (double)getdays(clim->year[rstep],clim->doy[rstep]);
				break;
			case 1:   // daily
				share->dayspan = 1;
				break;
			default:
				share->dayspan = (double)getdays(clim->year[rstep],clim->doy[rstep]);  //default for monthly
				break;
		}


		// Calculate cumulative GDD for other subroutines

		GDD = share->Tave * share->dayspan;

		share->TaveYr += GDD/365.0 ;
		share->PARYr += clim->par[rstep]* share->dayspan/365.0 ;

		if (GDD < 0 || clim->doy[rstep]<60 )  GDD = 0;   // need modification for tropical region

		share->GDDTot = share->GDDTot + GDD;


		// Phenology models altered by Aaron Teets (2021); aft49@nau.edu
		// Only for PnET-daily;
		// Four spring onset and two fall senescence models are included below,
		// new models are listed in increasing complexity;
		// Models were parameterized using PhenoCam data from 8 northern deciduous forests of North America


		// Select spring phenology subroutine here, models increase in complexity with selection values;

		//******************************************************

				spring_mod_select = 0;

		//******************************************************

		// For original PnET spring phenology model, use 0
		// For Temporal Time (TT) model, use 1
		// For M1 model, use 2
		// For PTT model, use 3
		// For Parallel model (PA), use 4
		// For Sequential M1 model (SM1), use 5
		// For Fixed spring phenology model, use 6


		// start phenology models

		if (spring_mod_select == 0)
		{

		// ORIGINAL PnET Phenology model

		F_crit1 = 100;
		F_crit2 = 900;

		if (clim->doy[rstep] < 60)  Tc = 0, share->Tsum = 0, share->CDDTot=0;
		Tc = share->Tave;

		}

		else if (spring_mod_select == 1)
		{
		// TEMPORAL TIME (TT) MODEL;
		// Temperature forcing only

		// parameters for TT model

				t0 = 92;
				T_base = 3.35;
				F_crit1 = 138.26;
				F_crit2 = 283.04;

				if (clim->doy[rstep] <= t0)  Tc = 0, share->Tsum = 0, share->CDDTot=0;
				else if(clim->doy[rstep] > t0)  Tc = share->Tave - T_base;
				if (Tc < 0) Tc = 0;
		}

		else if (spring_mod_select == 2)
		{
		//M1 MODEL;
		// Temperature forcing and photoperiod requirements

		//parameters for M1 model

				t0 = 92; //par[1]
				T_base = 6.87; //par[2]
				k = 4.79; //par[3]
				F_crit1 = 266.75; //par[4], threshold for bud burst
				F_crit2 = 812.94; // threshold for leaf completion

				if (clim->doy[rstep] <= t0)  Tc = 0, share->Tsum = 0, share->CDDTot=0;
				else if(clim->doy[rstep] > t0)  Tc = share->Tave - T_base;
				if (Tc < 0) Tc = 0;
				Tc = pow(((share->DayLength/3600)/10),k)*Tc;
		}

		else if (spring_mod_select == 3)
		{
		//PTT MODEL;
		// Temperature forcing and photoperiod requirements

		//parameters for PTT model

				t0 = 91; //par[1]
				T_base = 4.98; //par[2]
				F_crit1 = 55.94; //par[4], threshold for bud burst
				F_crit2 = 126.86; // threshold for leaf completion

				if (clim->doy[rstep] <= t0)  Tc = 0, share->Tsum = 0, share->CDDTot=0;
				else if(clim->doy[rstep] > t0)  Tc = share->Tave - T_base;
				if (Tc < 0) Tc = 0;
				Tc = ((share->DayLength/3600)/24)*Tc;
		}

		else if (spring_mod_select == 4)
		{
		// PARALLEL (PA) MODEL;
		// Temperature forcing, and chilling requirements
		// parameters for PA model

				t0 = 91;//par[1]
				t0_chill = 72;//par[2]
				T_base = 7.38;//par[3]
				T_opt = 4.24;//par[4]
				T_min = -0.87;//par[5]
				T_max = 10.72;//par[6]
				C_ini = 0.087;//par[7]
				F_crit1 = 9.84;//par[8]
				F_crit2 = 28.33;//par[8.5]
				C_req = 142.49;//par[9]

				//chilling triangular temperature response
				if(share->Tave <= T_min) Rc = 0;
				else if (share->Tave > T_min && share->Tave < T_opt) Rc = (share->Tave - T_min) / (T_opt - T_min);
				else if (share->Tave > T_opt && share->Tave < T_max) Rc = (share->Tave - T_max) / (T_opt - T_max);
				else if (share->Tave >= T_max) Rc = 0;

				if (clim->doy[rstep] <= t0_chill)  Rc = 0, share->Sc = 0;
				share->Sc = share->Sc + Rc;

				k = C_ini + share->Sc * (1 - C_ini)/C_req;
				if(share->Sc > C_req) int k = 1;

				//forcing
				Tc = share->Tave - T_base;
				if (Tc < 0)  Tc = 0;
				Tc = Tc * k;
				if (clim->doy[rstep] <= t0)  Tc = 0, share->Tsum = 0, share->CDDTot=0;
		}

		else if (spring_mod_select == 5)
		{
		// SEQUENTIAL M1 (SM1) MODEL;
		// Temperature forcing, photoperiod, and chilling requirements
		// parameters for SM1 model

				t0 = 92;//par[1]
				t0_chill = 84;//par[2]
				T_base = 5.74;//par[3]
				T_opt = 9.95;//par[4]
				T_min = -1.34;//par[5]
				T_max = 14.85;//par[6]
				F_crit1 = 81.7;//par[7]
				F_crit2 = 189.8;//par[7.5]
				C_req = 4.36;//par[8]

				//chilling triangular temperature response
				if(share->Tave <= T_min) Rc = 0;
				else if (share->Tave > T_min && share->Tave < T_opt) Rc = (share->Tave - T_min) / (T_opt - T_min);
				else if (share->Tave > T_opt && share->Tave < T_max) Rc = (share->Tave - T_max) / (T_opt - T_max);
				else if (share->Tave >= T_max) Rc = 0;

				if (clim->doy[rstep] <= t0_chill)  Rc = 0;
				share->Sc = Rc;

				if(share->Sc > C_req) int k = 0;
				else int k = 1;

				//forcing
				Tc = share->Tave - T_base;
				if (Tc < 0)  Tc = 0;
				Tc = pow(((share->DayLength/3600) / 24), k) * Tc;
				if (clim->doy[rstep] <= t0_chill)  Tc = 0, share->Tsum = 0, share->CDDTot=0;
		}


		else if (spring_mod_select == 6)
				{
				// Fixed Phenology model
						share->Tsum = clim->doy[rstep];
						F_crit1 = 135;
						F_crit2 = 147;
						Tc = 0;
						if (clim->doy[rstep] < 60)  share->Tsum = 0, share->CDDTot=0;
		}



		// BUDBURST

		share->Tsum = share->Tsum + Tc;

		//if ((share->GDDTot >= veg->GDDFolStart) && (clim->doy[rstep] < veg->SenescStart))  // original
		//if ((share->GDDTot >= F_crit1) && (clim->doy[rstep] < veg->SenescStart - 70))
		if((share->Tsum >= F_crit1) && (clim->doy[rstep] < 210))
		{
			OldFolMass = share->FolMass;  
			GDDFolEff = (share->Tsum - F_crit1) / (F_crit2 - F_crit1);

			if (GDDFolEff < 0)GDDFolEff = 0;
			if (GDDFolEff > 1)GDDFolEff = 1;

			delGDDFolEff = GDDFolEff - share->OldGDDFolEff;
			share->FolMass = share->FolMass + (share->BudC * delGDDFolEff) / veg->CFracBiomass;
			share->FolProdCMo = (share->FolMass - OldFolMass) * veg->CFracBiomass;
			share->FolGRespMo = share->FolProdCMo * veg->GRespFrac;
			share->OldGDDFolEff = GDDFolEff;
		}
		else
		{
			share->FolProdCMo = 0; 
			share->FolGRespMo = 0;
		}
	}
	else
	{
		// Select fall phenology models here

			//******************************************************

						fall_mod_select = 0;

			//******************************************************

			// For Original PnET fall phenology model, use 0
			// For Cooling Degree Day (CDD) model, use 1
			// For Cooling Degree Day Photoperiod (CDDP) model, use 2
			// For Fixed fall phenology model, use 3


			// FALL SENESCENCE

			// Start Fall models

				share->FolLitM = 0;

				// Original PnET Fall model
				if (fall_mod_select == 0)
					{
					share->FolLitM = 0;
							if ((share->PosCBalMass < share->FolMass) && (clim->doy[rstep] > 268))
							{


								if ((share->PosCBalMass) > (veg->FolMassMin))
								{
									FolMassNew = share->PosCBalMass;
								}
								else
								{
									FolMassNew = veg->FolMassMin;
								}


								if (FolMassNew == 0)
								{
									share->LAI = 0;
								}
								else if (FolMassNew < share->FolMass)
								{
									share->LAI = share->LAI * (FolMassNew / share->FolMass);
								}
								if (FolMassNew < share->FolMass)
								{
									share->FolLitM = share->FolMass - FolMassNew;
								}
								share->FolMass = FolMassNew;
							}
						}



				// COOLING DEGREE DAY (CDD) FALL SENESCENCE MODEL;
				// Accumulated cold temperature forcing
				else if (fall_mod_select == 1)
					{
							fall_t0 = 226; //
							fall_T_base = 25.08; //
							fall_F_crit1 = -384.8; // threshold for start of leaf senescence
							fall_F_crit2 = -784.22; // threshold for end of leaf senescence

							CDD = (share->Tave - fall_T_base) * share->dayspan;

							if (CDD > 0) CDD=0;

							if (clim->doy[rstep] <= fall_t0)  CDD = 0;

							share->CDDTot = share->CDDTot + CDD;

							share->FolLitM = 0;

						if ((share->CDDTot <= fall_F_crit1) && (clim->doy[rstep] > fall_t0))
						{
							if (share->CDDTot > fall_F_crit2)
							{
									FolMassNew = share->FolMass * (fall_F_crit2 - share->CDDTot) / (fall_F_crit2 - fall_F_crit1);
							}
						else if (share->CDDTot < fall_F_crit2)
						{
							FolMassNew = veg->FolMassMin;
						}
							if (FolMassNew == 0)
							{
								share->LAI = 0;
							}
							else if (FolMassNew < share->FolMass)
							{
								share->LAI = share->LAI * (FolMassNew / share->FolMass);
							}
						if (FolMassNew < share->FolMass)
						{
								share->FolLitM = share->FolMass - FolMassNew;
						}
								share->FolMass = FolMassNew;
						}
					}


				//  COOLING DEGREE DAY PHOTOPERIOD (CDDP) FALL SENESCENCE MODEL;
				// Accumulated cold temperature forcing
				else if (fall_mod_select == 2)
					{
							fall_t0 = 250; //
							fall_T_base = 10.68; //
							fall_F_crit1 = 27.12; // threshold for start of leaf senescence
							fall_F_crit2 = 96.46; // threshold for end of leaf senescence

								//fall_tmin = float clim->tmin;

								CDD = (clim->tmin[rstep] - fall_T_base) * share->dayspan;

								if (CDD > 0) CDD=0;
								CDD = ((1 - (share->DayLength/3600))/24) * CDD;
								if (clim->doy[rstep] <= fall_t0)  CDD = 0;

								share->CDDTot = share->CDDTot + CDD;

								share->FolLitM = 0;

								if ((share->CDDTot >= fall_F_crit1) && (clim->doy[rstep] > fall_t0))
								{
									if (share->CDDTot < fall_F_crit2)
									{
											FolMassNew = share->FolMass * (fall_F_crit2 - share->CDDTot) / (fall_F_crit2 - fall_F_crit1);
									}
								else if (share->CDDTot > fall_F_crit2)
								{
									FolMassNew = veg->FolMassMin;
								}
									if (FolMassNew == 0)
									{
										share->LAI = 0;
									}
									else if (FolMassNew < share->FolMass)
									{
										share->LAI = share->LAI * (FolMassNew / share->FolMass);
									}
								if (FolMassNew < share->FolMass)
								{
										share->FolLitM = share->FolMass - FolMassNew;
								}
										share->FolMass = FolMassNew;
								}
						}

				// Fixed Fall model
				else if (fall_mod_select == 3)
							{
							f_fixed1= 265; // DOY to start fall transition
							f_fixed2= 280; // DOY to complete fall transition

							if ((clim->doy[rstep] > f_fixed1))
								{
									if (clim->doy[rstep] <= f_fixed2)
									{
										FolMassNew = share->FolMass * ((f_fixed2 - float(clim->doy[rstep]))/(f_fixed2-f_fixed1));
									}
									else if (clim->doy[rstep] > f_fixed2)
									{
										FolMassNew = veg->FolMassMin;
									}
							if (FolMassNew == 0)
							{
								share->LAI = 0;
							}
							else if (FolMassNew < share->FolMass)
							{
								share->LAI = share->LAI * (FolMassNew / share->FolMass);
							}
							if (FolMassNew < share->FolMass)
							{
								share->FolLitM = share->FolMass - FolMassNew;
							}
							share->FolMass = FolMassNew;
								}
							}
			}
		}

