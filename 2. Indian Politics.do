//Creating Political Controls//

//Global Path laptop//
global politics "C:\Users\Admin\Dropbox\Saptashya's thesis\data\ias tcpd\ias-release-main\Politics"
global election "C:\Users\Admin\Dropbox\Saptashya's thesis\data\ias tcpd\ias-release-main\Election Data"


//global paths ***When used in LAB***//
global election "C:\Users\sg575\Dropbox\Saptashya's thesis\data\ias tcpd\ias-release-main\Election Data"
global politics "C:\Users\sg575\Dropbox\Saptashya's thesis\data\ias tcpd\ias-release-main\Politics"



//Using Besley Burgess Politics Data//
use "$politics\Politics.dta", clear 

**Cleaning the data set and changing the cadre codes**
drop if year<1951
drop state
rename statenm cadre
gen cadre_code=1 if cadre=="Andhra Pradesh"
replace cadre_code=4 if cadre=="Bihar"
replace cadre_code=5 if cadre=="Chhattisgarh"
replace cadre_code=7 if cadre=="Gujarat"
replace cadre_code=8 if cadre=="Haryana"
replace cadre_code=9 if cadre=="Himachal Pradesh"
replace cadre_code=11 if cadre=="Jharkhand"
replace cadre_code=10 if cadre=="Jammu & Kashmir"
replace cadre_code=12 if cadre=="Karnataka"
replace cadre_code=13 if cadre=="Kerala"
replace cadre_code=14 if cadre=="Madhya Pradesh"
replace cadre_code=15 if cadre=="Maharashtra"
replace cadre_code=20 if cadre=="Orissa"
replace cadre_code=21 if cadre=="Punjab"
replace cadre_code=22 if cadre=="Rajasthan"
replace cadre_code=24 if cadre=="Tamil Nadu"
replace cadre_code=25 if cadre=="Telangana"
replace cadre_code=27 if cadre=="Uttar Pradesh"
replace cadre_code=28 if cadre=="Uttarakhand"
replace cadre_code=29 if cadre=="West Bengal"


order cadre_code year

rename cadre_code splitcadrecode
//Making the Politics cleaned monthly//

expand = 12
sort splitcadrecode year
bysort splitcadrecode year : gen month = _n
order splitcadrecode year month
save "$politics\Politics cleaned.dta", replace


//Making the Political infor post 2002 monthly//
use "$politics\Political infor post 2002.dta", clear
rename month e_month
expand = 12
sort splitcadrecode year
bysort splitcadrecode year : gen month = _n
replace e_month=. if month != e_month
order cadre splitcadrecode year month e_month
drop e_month
save "$politics\Political infor post 2002 monthly.dta", replace


//Appending the Politics cleaned monthly & Political infor post 2002 monthly data set//
use "$politics\Politics cleaned.dta", clear
merge m:m splitcadrecode year using "$politics\Political infor post 2002 monthly.dta"
sort splitcadrecode year month
drop elecdum elec_dum
drop _merge
save "$politics\Political all monthly.dta", replace


//Merging Election Panel with politics cleaned to form political control Data set//
use "$election\Electoral cycle replication 2.dta", clear
merge 1:1 splitcadrecode year month using "$politics\Political all monthly.dta"
drop _merge
order cadre year month splitcadrecode elec_dum 


/*
drop if splitcadrecode==.
drop _merge
save "$politics\Political Control.dta"

clear

//Merging political control with politic info post 2002 to form political control cleaned Data set//
use "$politics\Political Control.dta", clear 
merge 1:1 year month splitcadrecode elec_dum turnout inc incu ics cpi cpm bjp jd jp lkdp psp sp tdp agp jknc shs uc sad dmk r oth ind convote presrule no_seat no_cand using "$politics\Political infor post 2002.dta"
sort splitcadrecode year
drop if year[_n]==year[_n-1] & _merge==1
drop if year[_n]==year[_n+1] & _merge==1

*/

**Creating a unique number for making the parties vote continious**
sort splitcadrecode year month
drop if cadre=="Sikkim" & prev_elec_months== prev_elec_months[_n-1] //Deleting the duplicates for sikkim
gen ui = year* e_month* splitcadrecode
order cadre year month splitcadrecode elec_dum e_month ui
duplicates tag ui, gen (dup)
tab dup
replace ui=ui+1 if splitcadrecode==7 & ((year==2012 & month==12) | (year ==2017 & month==12))
replace ui=ui+1 if splitcadrecode==8 & ((year==1957 & month==5) | (year ==2019 & month==10))
drop dup
duplicates tag ui, gen (dup)
tab dup
drop dup
by splitcadrecode: replace ui=ui[_n-1] if missing(ui)

**Making The vote continious**
replace inc=inc[_n-1] if ui[_n-1]==ui[_n]
replace cpi = cpi[_n-1] if ui[_n-1]==ui[_n]
replace cpm = cpm[_n-1] if ui[_n-1]==ui[_n]
replace bjp = bjp[_n-1] if ui[_n-1]==ui[_n]
replace jd = jd[_n-1] if ui[_n-1]==ui[_n]
replace jp = jp[_n-1] if ui[_n-1]==ui[_n]
replace tdp = tdp[_n-1] if ui[_n-1]==ui[_n]
replace jknc = jknc[_n-1] if ui[_n-1]==ui[_n]
replace shs = shs[_n-1] if ui[_n-1]==ui[_n]
replace sad = sad[_n-1] if ui[_n-1]==ui[_n]
replace dmk = dmk[_n-1] if ui[_n-1]==ui[_n]
replace oth = oth[_n-1] if ui[_n-1]==ui[_n]
replace ind = ind[_n-1] if ui[_n-1]==ui[_n]
replace ncp = ncp[_n-1] if ui[_n-1]==ui[_n]
replace bjd = bjd[_n-1] if ui[_n-1]==ui[_n]
replace pdp = pdp[_n-1] if ui[_n-1]==ui[_n]
replace nc = nc[_n-1] if ui[_n-1]==ui[_n]
replace tmc = tmc[_n-1] if ui[_n-1]==ui[_n]
replace jmm = jmm[_n-1] if ui[_n-1]==ui[_n]
replace abrsp = abrsp[_n-1] if ui[_n-1]==ui[_n]
replace sdf = sdf[_n-1] if ui[_n-1]==ui[_n]
replace sjp = sjp[_n-1] if ui[_n-1]==ui[_n]
replace sc = sc[_n-1] if ui[_n-1]==ui[_n]
replace skc = skc[_n-1] if ui[_n-1]==ui[_n]
replace smp = smp[_n-1] if ui[_n-1]==ui[_n]
replace bsp = bsp[_n-1] if ui[_n-1]==ui[_n]
replace ljp = ljp[_n-1] if ui[_n-1]==ui[_n]
replace rdj = rdj[_n-1] if ui[_n-1]==ui[_n]
replace trs = trs[_n-1] if ui[_n-1]==ui[_n]
replace AIMIM = AIMIM[_n-1] if ui[_n-1]==ui[_n]
replace TRS = TRS[_n-1] if ui[_n-1]==ui[_n]
replace YSRC = YSRC[_n-1] if ui[_n-1]==ui[_n]
replace inld = inld[_n-1] if ui[_n-1]==ui[_n]
replace aiadmk = aiadmk[_n-1] if ui[_n-1]==ui[_n]





*Re indexing few party names**
rename  trs brs
rename rdj rjd
rename TRS trs
rename YSRC ysrc

**Replacing Missing values by 0**
replace inc=0 if missing(inc)
replace cpi=0 if missing(cpi)
replace cpm=0 if missing(cpm)
replace bjp=0 if missing(bjp)
replace jd=0 if missing(jd)
replace jp=0 if missing(jp)
replace tdp=0 if missing(tdp)
replace jknc=0 if missing(jknc)
replace shs=0 if missing(shs)
replace sad=0 if missing(sad)
replace dmk=0 if missing(dmk)
replace oth=0 if missing(oth)
replace ind=0 if missing(ind)
replace ncp=0 if missing(ncp)
replace bjd=0 if missing(bjd)
replace pdp=0 if missing(pdp)
replace nc=0 if missing(nc)
replace tmc=0 if missing(tmc)
replace jmm=0 if missing(jmm)
replace abrsp=0 if missing(abrsp)
replace sdf=0 if missing(sdf)
replace sjp=0 if missing(sjp)
replace sc=0 if missing(sc)
replace smp=0 if missing(smp)
replace bsp=0 if missing(bsp)
replace ljp=0 if missing(ljp)
replace rjd=0 if missing(rjd)
replace brs=0 if missing(brs)
replace inc=0 if missing(inc)
replace AIMIM=0 if missing(AIMIM)
replace trs=0 if missing(trs)
replace ysrc=0 if missing(ysrc)
replace inld=0 if missing(inld)
replace aiadmk=0 if missing(aiadmk)
replace incu=0 if missing(incu)
replace ics=0 if missing(ics)
replace lkdp=0 if missing(lkdp)
replace agp=0 if missing(agp)
replace psp=0 if missing(psp)
replace sp=0 if missing(sp)
replace uc=0 if missing(uc)
replace r=0 if missing(r)
drop skc

**Labeling new Parties**
label variable ncp "No_of_seats, Nationalist congress party"
label variable bjd "No_of_seats, Biju Janta dal"
label variable pdp "No_of_seats, J&K Peoples democratic party"
label variable nc "No_of_seats, J&K National Conference"
label variable tmc "No_of_seats, Trinamool Congress"
label variable jmm "No_of_seats, Jharkhand mukti morcha"
label variable abrsp "No_of_seats, Akhil Bharatiya Ram Rajya Parishad"
label variable sdf "No_of_seats, Sikkim Democratic Front"
label variable sjp "No_of_seats, Sikkim Sangram Parishad"
label variable sc "No_of_seats, Sikkim Congress"
label variable smp "No_of_seats, Samajwadi Party"
label variable bsp "No_of_seats, Bahujan Samaj Party"
label variable ljp "No_of_seats, Lok Janashakti Party"
label variable rjd "No_of_seats, Rashtriya Janta Dal"
label variable brs "No_of_seats, Bharatiya Rashtriya Samiti"
label variable AIMIM "No_of_seats, All India Majlis-e-Ittehadul Muslimeen"
label variable trs "No_of_seats, Telengana Rashtriya Samiti"
label variable ysrc "No_of_seats, Yuvajana Shramika Rythu Congress Party"
label variable inld "No_of_seats, Indian National Lok Dal"
label variable aiadmk "No_of_seats, All India Anna Dravida Munnetra Kazhagam"

**Generating total seats of tbe parties**
gen inc_tot= inc
gen left_tot= ics + cpi + cpm + psp + sp
gen bjp_tot= bjp + jp + abrsp
gen oth_tot= incu + bjd + jd + lkdp + tdp + agp + jknc + shs + uc + sad + dmk + r + oth + ind + ncp + pdp + nc + tmc + jmm + sdf + sjp + sc + smp + bsp + ljp + rjd + brs + AIMIM + trs + ysrc + inld + aiadmk
gen overall_tot = inc_tot + left_tot + bjp_tot + oth_tot

**Generating Vote shares**
gen inc_share = inc_tot/overall_tot
gen left_share = left_tot/overall_tot
gen bjp_share = bjp_tot/overall_tot
gen oth_share = oth_tot/overall_tot

save "$politics\Political Control.dta", replace

clear




use "$politics\Political Control.dta", clear
keep cadre year month splitcadrecode elec_dum e_month prev_elec_months next_elec_months inc_share left_share bjp_share oth_share
order cadre splitcadrecode year
label variable inc_share "Indian National Congress Vote Share"
label variable left_share "Left Parties Vote Share"
label variable bjp_share "BJP Right Wings and Hindutva parties Vote Share"
label variable oth_share "other regional parties Vote Share"

save "$politics\Political Control Cleaned.dta", replace




















