/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "surfaceChemistryModel.H"
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::surfaceChemistryModel
(
    ReactionThermo& thermo
)
:
    BasicChemistryModel<ReactionThermo>(thermo),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),
    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_
    (
        BasicChemistryModel<ReactionThermo>::template getOrDefault<scalar>
        (
            "Treact",
            0.0
        )
    ),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_)
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
            )
        );
    }

    Info<< "homogenous chemistry : Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
    Info<< "heterogenous chemistry : Number of species = " << nSpecieSurf_
	<< " and surface reaction s = " << nReactionSurf_ << endl;
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::
~surfaceChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    dcdt = Zero;

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            dcdt[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            dcdt[si] += sr*omegai;
        }
    }
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const Reaction<ThermoType>& R = reactions_[index];
    scalar w = omega(R, c, T, p, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c);
    const scalar kr = R.kr(kf, p, T, c);

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0.0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0.0), exp);
        }
    }
    cf = max(c[lRef], 0.0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0.0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(c[si], 0.0), exp);
        }
    }
    cr = max(c[rRef], 0.0);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class ReactionThermo, class ThermoType>
void Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0.0);
    }

    omega(c_, T, p, dcdt);

    // Constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar W = specieThermo_[i].W();
        rho += W*c_[i];
    }
    scalar cp = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        cp += c_[i]*specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    scalar dT = 0.0;
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar hi = specieThermo_[i].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class ReactionThermo, class ThermoType>
void Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0.0);
    }

    dfdc = Zero;

    // Length of the first argument must be nSpecie_
    omega(c_, T, p, dcdt);

    forAll(reactions_, ri)
    {
        const Reaction<ThermoType>& R = reactions_[ri];

        const scalar kf0 = R.kf(p, T, c_);
        const scalar kr0 = R.kr(kf0, p, T, c_);

        forAll(R.lhs(), j)
        {
            const label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c_[si] > SMALL)
                        {
                            kf *= el*pow(c_[si], el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c_[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c_[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc(si, sj) -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc(si, sj) += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            const label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar er = R.rhs()[i].exponent;
                if (i == j)
                {
                    if (er < 1.0)
                    {
                        if (c_[si] > SMALL)
                        {
                            kr *= er*pow(c_[si], er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c_[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= pow(c_[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc(si, sj) += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc(si, sj) -= sr*kr;
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1.0e-3;

    omega(c_, T + delta, p, dcdt_);
    for (label i=0; i<nSpecie_; i++)
    {
        dfdc(i, nSpecie_) = dcdt_[i];
    }

    omega(c_, T - delta, p, dcdt_);
    for (label i=0; i<nSpecie_; i++)
    {
        dfdc(i, nSpecie_) = 0.5*(dfdc(i, nSpecie_) - dcdt_[i])/delta;
    }

    dfdc(nSpecie_, nSpecie_) = 0;
    dfdc(nSpecie_ + 1, nSpecie_) = 0;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("small", dimTime, SMALL),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
   
    // ajout : calcul le temps de reaction chimique pour 
    // des reactions en volume + surface

    volScalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const volScalarField& rho = trho();
    const volScalarField& T   = this->thermo().T();
    const volScalarField& p   = this->thermo().p();
    const label nReaction     = reactions_.size();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    if (this->chemistryVol_)
    {
    	forAll(rho, celli)
        {
        	const scalar rhoi = rho[celli];
            	const scalar Ti = T[celli];
            	const scalar pi = p[celli];

            	scalar cSum = 0.0;

            	for (label i=0; i<nSpecie_; i++)
            	{
                	c_[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
                	cSum += c_[i];
            	}

            	forAll(reactions_, i)
            	{
                	const Reaction<ThermoType>& R = reactions_[i];

                	omega(R, c_, Ti, pi, pf, cf, lRef, pr, cr, rRef);

                	forAll(R.rhs(), s)
                	{
                    		tc[celli] += R.rhs()[s].stoichCoeff*pf*cf;
                	}
            	}
            	tc[celli] = nReaction*cSum/tc[celli];
        }
    }
    else if (this->chemistrySurf_)
    {
	label patchID = this->mesh().boundaryMesh().findPatchID("catalyticWall");
    	fvPatchScalarField  rhoPatchI  =  rho.boundaryField()  [patchID];
	fvPatchScalarField  TPatchI    =  T.boundaryField()    [patchID];
	fvPatchScalarField  pPatchI    =  p.boundaryField()    [patchID];
	fvPatchScalarField&  tcPatchI  =  tc.boundaryFieldRef()[patchID];
        forAll(rhoPatchI, faceI)
        {
      		const scalar rhoi = rhoPatchI[faceI];
            	const scalar Ti   = TPatchI  [faceI];
            	const scalar pi   = pPatchI  [faceI];

            	scalar cSum = 0.0;

            	for (label i=0; i<nSpecie_; i++)
            	{
			fvPatchScalarField Yi = Y_[i].boundaryField()[patchID];    
                	c_[i] = rhoi*Yi[faceI]/specieThermo_[i].W();
                	cSum += c_[i];
            	}

            	forAll(reactions_, i)
            	{
               		const Reaction<ThermoType>& R = reactions_[i];
                	omega(R, c_, Ti, pi, pf, cf, lRef, pr, cr, rRef);

                	forAll(R.rhs(), s)
                	{
                   		tcPatchI[faceI] += R.rhs()[s].stoichCoeff*pf*cf;
                	}
            	}
            	tcPatchI[faceI] = nReaction*cSum/tc[faceI];
        }
    }

    // ttc.ref().correctBoundaryConditions();
  
    return ttc;
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::Qdot
(
) const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    volScalarField& Qdot = tQdot.ref();

    // ----------------- AJOUT : REACTION DE SURFACE --------- // 
    if (this->chemistryVol_)
    {
        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermo_[i].Hc();
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }
    else if (this->chemistrySurf_)
    {
	    	label patchID = this->mesh().boundaryMesh().findPatchID("catalyticWall");
		fvPatchScalarField& QdotPatchI = Qdot.boundaryFieldRef()[patchID];
		forAll(Y_,i)
		{

                	const scalar hi    = specieThermo_[i].Hc();
			fvPatchScalarField RRPatchI    = RR_[i].boundaryField()[patchID];
			fvPatchScalarField YPatchI     = Y_[i].boundaryField()[patchID];
			forAll(QdotPatchI,faceI)
			{	
				const scalar x = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].x();
				const scalar y = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].y();

				if(mag(x) <= this->catalyticLength_ && mag(y)<= this->catalyticLength_)
				{
                			QdotPatchI[faceI] -= hi*RRPatchI[faceI];
				}
				else
				{
					QdotPatchI[faceI] = 0.;
				}

			}
		}	
    }

    return tQdot;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );
 
    volScalarField& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const volScalarField& rho = trho();

    const volScalarField& T = this->thermo().T();
    const volScalarField& p = this->thermo().p();

    if(this->chemistryVol_)
    {
    	forAll(rho, celli)
    	{
       		const scalar rhoi = rho[celli];
        	const scalar Ti = T[celli];
        	const scalar pi = p[celli];
	
        	for (label i=0; i<nSpecie_; i++)
        	{
        	    const scalar Yi = Y_[i][celli];
        	    c_[i] = rhoi*Yi/specieThermo_[i].W();
        	}

        	const scalar w = omegaI
        	(
            		ri,
            		c_,
            		Ti,
            		pi,
            		pf,
            		cf,
           		lRef,
            		pr,
            		cr,
            		rRef
        	);

        	RR[celli] = w*specieThermo_[si].W();
	}
    }
    //-------------- AJOUT ----------- //
    else if(this->chemistrySurf_)
    {
	label patchID = this->mesh().boundaryMesh().findPatchID("catalyticWall");
	fvPatchScalarField  rhoPatchI  = rho.boundaryField()  [patchID];
	fvPatchScalarField  TPatchI    = T.boundaryField()    [patchID];
	fvPatchScalarField  pPatchI    = p.boundaryField()    [patchID];
	fvPatchScalarField&  RRPatchI  = RR.boundaryFieldRef()[patchID];
	forAll(rhoPatchI, faceI)
    	{
       		const scalar rhoi = rhoPatchI[faceI];
        	const scalar Ti   = TPatchI[faceI];
        	const scalar pi   = pPatchI[faceI];
	
		const scalar x = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].x();
		const scalar y = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].y();
        	for (label i=0; i<nSpecie_; i++)
        	{
        	    const fvPatchScalarField Yi = Y_[i].boundaryField()[patchID];
        	    c_[i] = rhoi*Yi[faceI]/specieThermo_[i].W();
        	}

        	const scalar w = omegaI
        	(
            		ri,
            		c_,
            		Ti,
            		pi,
            		pf,
            		cf,
           		lRef,
            		pr,
            		cr,
            		rRef
        	);

		if(mag(x) <= this->catalyticLength_ && mag(y) <= this->catalyticLength_)
		{
        		RRPatchI[faceI] = w*specieThermo_[si].W();
		}
		else
		{
			RRPatchI[faceI] = 0.;
		}
	}
    }
    return tRR;
}

template<class ReactionThermo, class ThermoType>
void Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::calculate()
{
    if (!this->chemistryVol_ && !this->chemistrySurf_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const volScalarField& rho = trho();
    const volScalarField& T   = this->thermo().T();
    const volScalarField& p   = this->thermo().p();

    if(this->chemistryVol_)
    {
  	forAll(rho, celli)
    	{
       		const scalar rhoi = rho[celli];
        	const scalar Ti = T[celli];
       		const scalar pi = p[celli];
        	for (label i=0; i<nSpecie_; i++)
        	{
            		const scalar Yi = Y_[i][celli];
            		c_[i] = rhoi*Yi/specieThermo_[i].W();
        	}

        	omega(c_, Ti, pi, dcdt_);

        	for (label i=0; i<nSpecie_; i++)
        	{
            		RR_[i][celli] = dcdt_[i]*specieThermo_[i].W();
        	}
    	}
    }
    //------------- AJOUT : REACTION DE SURFACE SUR UNE BC --------------//
    else if(this->chemistrySurf_)
    {
	label patchID = this->mesh().boundaryMesh().findPatchID("catalyticWall");
	fvPatchScalarField TPatchI    = T.boundaryField()  [patchID];
	fvPatchScalarField pPatchI    = p.boundaryField()  [patchID];
	fvPatchScalarField rhoPatchI  = rho.boundaryField()[patchID];
	forAll(rhoPatchI, faceI)
    	{
       		const scalar rhoi = rhoPatchI[faceI];
        	const scalar Ti   = TPatchI  [faceI];
       		const scalar pi   = pPatchI  [faceI];
	
		const scalar x = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].x();
		const scalar y = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].y();

		// calcul pour chaque face la concentration ci
        	for (label i=0; i<nSpecie_; i++)
        	{
			fvPatchScalarField YPatchI = Y_[i].boundaryField()[patchID];
            		const scalar Yi            = YPatchI[faceI];
            		c_[i] = rhoi*Yi/specieThermo_[i].W();
        	}

        	omega(c_, Ti, pi, dcdt_);

        	for (label i=0; i<nSpecie_; i++)
        	{
			fvPatchScalarField& RRPatchI = RR_[i].boundaryFieldRef()
				[patchID];
			if(mag(x) <= this->catalyticLength_ && 
					mag(y) <= this->catalyticLength_)
			{
            			RRPatchI[faceI] = dcdt_[i]*specieThermo_[i].W();
			}
			else
			{
            			RRPatchI[faceI] = 0.;
			}
        	}
    	}    

    }
}

template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;

    if ( (!this->chemistryVol_) && (!this->chemistrySurf_) )
    {
        return deltaTMin;
    }

    // modification : considere les champs de T,p,rho comme
    // des volScalarField -> permet d'acceder aux CL
 
    tmp<volScalarField> trho(this->thermo().rho());
    const volScalarField& rho = trho();

    const volScalarField& T = this->thermo().T();
    const volScalarField& p = this->thermo().p();

    scalarField c0(nSpecie_);

    if(this->chemistryVol_)
    {
    	forAll(rho, celli)
    	{
       		scalar Ti = T[celli];

        	if (Ti > Treact_)
        	{
            		const scalar rhoi = rho[celli];
            		scalar pi = p[celli];

            		for (label i=0; i<nSpecie_; i++)
            		{
                		c_[i] = rhoi*Y_[i][celli]/
					specieThermo_[i].W();
                		c0[i] = c_[i];
            		}

            		// Initialise time progress
            		scalar timeLeft = deltaT[celli];

            		// Calculate the chemical source terms
            		while (timeLeft > SMALL)
            		{
                		scalar dt = timeLeft;
                		this->solve(c_, Ti, pi, dt, 
						this->deltaTChem_[celli]);
                		timeLeft -= dt;
            		}

            		deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            		this->deltaTChem_[celli] =
                	min(this->deltaTChem_[celli], this->deltaTChemMax_);

            		for (label i=0; i<nSpecie_; i++)
            		{
                		RR_[i][celli] =(c_[i] - c0[i])*
					specieThermo_[i].W()/deltaT[celli];
            		}
        	}
        	else
        	{
            		for (label i=0; i<nSpecie_; i++)
            		{
                		RR_[i][celli] = 0;
            		}
        	}
    	}
    }

    // ------------------- AJOUT : REACTION A LA SUFACE SUR UNE BC --------//
    else if( this->chemistrySurf_)
    {
	const label patchID                  = this->mesh().
		boundaryMesh().findPatchID("catalyticWall");
	const fvPatchScalarField TPatchI     = T.boundaryField()  [patchID];
	const fvPatchScalarField rhoPatchI   = rho.boundaryField()[patchID];
	const fvPatchScalarField pPatchI     = p.boundaryField()  [patchID];
	fvPatchScalarField& deltaTChemPatchI = this->deltaTChem_.
		boundaryFieldRef()[patchID];
	forAll(rhoPatchI, faceI)
    	{
		const scalar x = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].x();
		const scalar y = this->mesh().boundaryMesh()[patchID].faceCentres()[faceI].y();
       		scalar Ti = TPatchI[faceI];
        	if (Ti > Treact_)
        	{
            		const scalar rhoi = rhoPatchI[faceI];
            		scalar pi         = pPatchI[faceI];
            		for (label i=0; i<nSpecie_; i++)
            		{
				fvPatchScalarField& Yi = Y_[i].boundaryFieldRef()[patchID];
				if(Yi[faceI]<0.)
				{
					Yi[faceI] = 0.;
				}
                		c_[i] = rhoi*Yi[faceI]/specieThermo_[i].W();
                		c0[i] = c_[i];
            		}

            		// Initialise time progress
            		scalar timeLeft = deltaT[faceI];

            		// Calculate the chemical source terms
            		while (timeLeft > SMALL)
            		{
                		scalar dt = timeLeft;
                		this->solve(c_, Ti, pi, dt,deltaTChemPatchI[faceI]);
                		timeLeft -= dt;
            		}

            		deltaTChemPatchI[faceI] =
                	min(deltaTChemPatchI[faceI], this->deltaTChemMax_);

            		for (label i=0; i<nSpecie_; i++)
            		{
				fvPatchScalarField& RRiPatchI = RR_[i].
					boundaryFieldRef()[patchID];
				if(mag(x) <= this->catalyticLength_ && mag(y) <= this->catalyticLength_)	
				{		
                			RRiPatchI[faceI] =(c_[i] - c0[i])*specieThermo_[i].W()
					/deltaT[faceI];
				}
				else
				{
					RRiPatchI[faceI] = 0.;
				}
            		}
        	}
        	else
        	{
            		for (label i=0; i<nSpecie_; i++)
            		{
				fvPatchScalarField& RRiPatchI = RR_[i].
					boundaryFieldRef()[patchID];
                		RRiPatchI[faceI] = 0;
            		}
		}
    	}
    }
    return deltaTMin;
}

template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}

template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::surfaceChemistryModel<ReactionThermo, ThermoType>::solve
(
    const volScalarField& deltaT
)
{
    return this->solve<volScalarField>(deltaT);
}


// ************************************************************************* //
