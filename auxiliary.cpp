#include "auxiliary.h"

string recordFile(string path)
{
	time_t t = time(NULL);
	struct tm current_time;
	localtime_s(&current_time, &t);
	string directory = path + to_string(current_time.tm_year + 1900) + string("-")
		+ to_string(current_time.tm_mon + 1) + string("-")
		+ to_string(current_time.tm_mday) + string("(")
		+ to_string(current_time.tm_hour) + string("-")
		+ to_string(current_time.tm_min) + string(")");

	string mkdir = string("md ") + directory;
	if (system(mkdir.c_str()))
	{
		cerr << "DATA DIRECTORY CREATING ERROR! " << endl;
		system("pause");
	}
	return directory;
}

ofstream logFile(string directory)
{
	string logfile = directory + string("\\simulation_log.txt");
	ofstream log_file(logfile.c_str());
	if (!log_file.is_open())
	{
		cerr << " LOG FILE OPEN ERROR !" << endl;
		system("pause");
	}
	log_file << directory << endl;
	return log_file;
}

void recordInitial(Species& electron, Species& ion, ofstream& log_file)
{
	log_file << "initial electron kinetic energy: " << showpoint << electron.kinetic() << endl;
	log_file << "initial ion kinetic energy: " << showpoint << ion.kinetic() << endl;
	log_file << "initial electron number: " << electron.part.size() << endl;
	log_file << "initial ion number: " << ion.part.size() << endl;
	cout << "initial electron kinetic energy: " << showpoint << electron.kinetic() << endl;
	cout << "initial ion kinetic energy: " << showpoint << ion.kinetic() << endl;
	cout << "initial electron number: " << electron.part.size() << endl;
	cout << "initial ion number: " << ion.part.size() << endl;
	log_file << endl;
	cout << endl;
}

void recordLog(double INITIAL_PARTICLE, int main_iter, Species& electron, Species& ion, double totalElectKinetic, double totalIonKinetic, ofstream& log_file)
{
	/* check electron */
	log_file << "Main loop number:" << main_iter << ",electron kinetic energy " << 100.0 * electron.kinetic() / totalElectKinetic << " %" << endl;
	log_file << "Main loop number:" << main_iter << ",electron number contains " << 100.0 * electron.part.size() / INITIAL_PARTICLE << " %" << endl;
	cout << "Main loop number:" << main_iter << ",electron kinetic energy " << 100.0 * electron.kinetic() / totalElectKinetic << " %" << endl;
	cout << "Main loop number:" << main_iter << ",electron number contains " << 100.0 * electron.part.size() / INITIAL_PARTICLE << " %" << endl;

	/* check ion */
	log_file << "Main loop number:" << main_iter << ",ion kinetic energy " << 100.0 * ion.kinetic() / totalIonKinetic << " %" << endl;
	log_file << "Main loop number:" << main_iter << ",ion number contains " << 100.0 * ion.part.size() / INITIAL_PARTICLE << " %" << endl;
	cout << "Main loop number:" << main_iter << ",ion kinetic energy " << 100.0 * ion.kinetic() / totalIonKinetic << " %" << endl;
	cout << "Main loop number:" << main_iter << ",ion number contains " << 100.0 * ion.part.size() / INITIAL_PARTICLE << " %" << endl;
	log_file <<"initialization complete!" << endl;
	cout <<"initialization complete!" << endl;
}

void recordData(int main_iter, string directory, Species& electron, Species& ion, vectorField& electCurrent,
	scalarField& electDen, vectorField& ionCurrent, scalarField& ionDen, scalarField& phi, scalarField& rho, vectorField& E_l,
	vectorField& E_t,vectorField& B )
{
	string path = directory + string("\\") + to_string(main_iter);

	/* save field data */
	electCurrent.saveInfo(path + string("Jex.dat"),
		                  path + string("Jey.dat"),
		                  path + string("Jez.dat"));
	electDen.saveInfo(path + string("electDen.dat"));
	ionCurrent.saveInfo(path + string("Jix.dat"),
		                path + string("Jiy.dat"),
		                path + string("Jiz.dat"));
	ionDen.saveInfo(path + string("ionDen.dat"));
	phi.saveInfo(path + string("phi.dat"));
	rho.saveInfo(path + string("rho.dat"));
	E_l.saveInfo(path + string("E_lx.dat"), 
		         path + string("E_ly.dat"),
		         path + string("E_lz.dat"));
	E_t.saveInfo(path + string("E_tx.dat"),
		         path + string("E_ty.dat"),
		         path + string("E_tz.dat"));
	B.saveInfo(path + string("Bx.dat"),                                                                              
		       path + string("By.dat"),
		       path + string("Bz.dat"));

	/* save particle data */
	electron.saveParticle(path + string("electron.dat"));
	ion.saveParticle(path + string("ion.dat"));
	electron.saveGlobal(path + string("electron_global.dat"));
	ion.saveGlobal(path + string("ion_global.dat"));
}