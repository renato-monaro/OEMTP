/**
 *  Copyright 2015 by Renato Monaro and Silvio Giuseppe
 *
 * This file is part of Open Electromagnetic Transient Program - OEMTP.
 * 
 * OEMTP is free software: you can redistribute 
 * it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation, either 
 * version 3 of the License, or (at your option) any later version.
 * 
 * Some open source application is distributed in the hope that it will 
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @license GPL-3.0+ <http://spdx.org/licenses/GPL-3.0+>
 */

#ifndef INTERFACE_ORELAY_H
#define INTERFACE_ORELAY_H



#ifdef OPENRELAY
#include<relay.h>
#include "components.h"
enum {TYPE_CURRENT=0,TYPE_VOLTAGE,TYPE_VALUE,TYPE_DIGITAL};
namespace orelay{
class OEMTPAcquisition: public Acquisition{
	public:
		OEMTPAcquisition(double dT);

		void Join(Channel<analogic>*,oemtp::Component*,unsigned Type, unsigned Pos);
		void Join(Channel<analogic>*,oemtp::Component*,unsigned Type);
		void Join(Channel<analogic>*,oemtp::Ammeter*);
		void Join(Channel<analogic>*,oemtp::Voltmeter*);

		void Join(Channel<digital>*,oemtp::Switch*,unsigned Pos, bool direction); 
		void Join(Channel<digital>*,oemtp::Switch*,bool direction); 
		void Join(Channel<digital>*,oemtp::Switch*); 

		virtual bool Run();
		virtual bool Prepare(float); 
		void RefreshTime();
	protected:
		vector <oemtp::Component* > ComponentL;
		vector <oemtp::Switch* > SwitchL;
		vector <Channel<analogic>* > AnalogicRes;
		vector <Channel<digital>* > DigitalRes;
		vector <Channel<analogic>* > AnalogicSw;
		vector <Channel<digital>* > DigitalSw;
		vector <unsigned > TypeRes;
		vector <bool > DirectionRes;
		vector <unsigned > TypeSw;
		vector <bool > DirectionSw;
		vector <unsigned> Position;
		timer dTns;
		timer T;
	private:
		analogic Value;
	};
}
#endif
#endif /* !INTERFACE_ORELAY_H */
