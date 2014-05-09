var fs = require('fs');
/*
var data = new Array(70);

for (var i = 0; i < data.length; i++) {
	var folder = "./set" + i + "/";

	var motif_txt      = fs.readFileSync( folder + "motif.txt").toString();
	motif_txt          = motif_txt.match(/[^MOTIF](A|C|G|T|[*])(A|C|G|T|[*])*\/)[0].trim(); // change \/ to just /
	var runtime        = parseInt( fs.readFileSync( folder + "runtime.txt" ).toString().trim().replace("\n", "") );

	var predictedSites = fs.readFileSync( folder + "predictedsites.txt" ).toString().split("\n");
	predictedSites.pop();

	var sites = fs.readFileSync(folder + "sites.txt" ).toString().split("\n");
	sites.pop();
	for (var a = 0; a < sites.length; a++) {
		sites[a] = parseInt(sites[a]);
	};

	for (var a = 0; a < predictedSites.length; a++) {
		predictedSites[a] = parseInt(predictedSites[a]);
	};

	var predictedMotifs = fs.readFileSync( folder + "predictedmotif.txt" ).toString().split("\n");
	predictedMotifs.splice(0,1);
	predictedMotifs.splice(predictedMotifs.length-1,1);
	predictedMotifs.splice(predictedMotifs.length-1,1);



	data[i] = {
		setNumber       : i,
		correctMotif    : motif_txt,
		runtime         : runtime,
		predictedSites  : predictedSites,
		sites           : sites,
		predictedMotifs : predictedMotifs
	};
};

fs.writeFileSync("data.json", JSON.stringify(data, null, '\t') );
*/

var data_json = require('./data.json');
var program_output = require('./program_output.json');

function countCharacter(string, char) {
	var count = 0;
	for (var i = 0; i < string.length; i++) {
		count += char == string.charAt(i);
	};
	return count;
}

function isInArray (key, array) {
	for (var i = 0; i < array.length; i++) {
		if ( key == array[i] )
			return true;
	};
	return false;
}

var data = new Array(70);

for (var i = 0; i < data_json.length; i++) {
	data[i] = {
		MotifLength : data_json[i].correctMotif.length,
		Accuracy : data_json[i].predictedSitesPercentage,
		RelativeEntropy : program_output[i].relativeEntropy
	};
};

fs.writeFileSync("scatterplot_accuracy_vs_re.json", JSON.stringify(data, null, '\t') );
