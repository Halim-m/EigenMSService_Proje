package service;

import com.sun.tools.sjavac.Source;
import domain.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class EigenMSService extends AbstractService{
    @Override
    public String getServiceName() {
        return "EigenMSService";
    }

    @Override
    public DNMSO run(DNMSO targetDNMSO, String[] args) throws IOException {
        processSettings(targetDNMSO, args);
        if (getProperties().get("command").equals("read")) return read();
        return null;
    }

    private DNMSO read() {
        String predictionFilePath = getProperties().get(ServiceTag.PREDICTION_FILE_PATH.toString());
        File predictionFile = new File(predictionFilePath);
        DNMSO targetDNMSO = getDNMSO();
        String size = getProperties().get(ServiceTag.NUMBER_OF_PREDICTION.toString());

        Pattern pattern = Pattern.compile("^(\\S+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)");

        String fileContent;

        if (targetDNMSO == null) {
            DNMSOFactory dnmsoFactory = new DNMSOFactory();
            targetDNMSO = dnmsoFactory.createDNMSO();
        }

        String sampleName;
        Score scoreAbundanceAverage;
        Score scoreAbundanceStDev;
        Score scoreSLiCScoreAverage;
        Score scoreHighMSMSScore;
        Score scoreHighDiscriminantScore;
        Score scorePMTQualityScore;
        Collection<Spectrum> spectra = new LinkedList<>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(predictionFile));
            while ((fileContent = in.readLine()) != null && check(targetDNMSO.getPrediction().size(), size)) {
                Matcher matcher = pattern.matcher(fileContent);
                while (matcher.find() && check(targetDNMSO.getPrediction().size(), size)) {
                    Prediction prediction = new Prediction();
                    prediction.setScore(new LinkedList<>());
                    if(matcher.group(1) != null) {
                        sampleName = matcher.group(1);
                        Spectrum spectrum = new Spectrum();
                        spectrum.setSpectrumId(sampleName);
                        spectra.add(spectrum);
                    }
                    if(matcher.group(5) != null) {
                        Double abundanceAverage =Double.parseDouble(matcher.group(5));
                        scoreAbundanceAverage = new Score();
                        scoreAbundanceAverage.setScoreName("abundanceAverage");
                        scoreAbundanceAverage.setScoreValue(abundanceAverage);
                        prediction.getScore().add(scoreAbundanceAverage);
                    }
                    if(matcher.group(6) != null) {
                        Double abundanceStDev = Double.parseDouble(matcher.group(6));
                        scoreAbundanceStDev = new Score();
                        scoreAbundanceStDev.setScoreName("abundanceStDev");
                        scoreAbundanceStDev.setScoreValue(abundanceStDev);
                        prediction.getScore().add(scoreAbundanceStDev);
                    }
                    if(matcher.group(7) != null) {
                        Double slicScoreAverage = Double.parseDouble(matcher.group(7));
                        scoreSLiCScoreAverage = new Score();
                        scoreSLiCScoreAverage.setScoreName("slicScoreAverage");
                        scoreSLiCScoreAverage.setScoreValue(slicScoreAverage);
                        prediction.getScore().add(scoreSLiCScoreAverage);
                    }
                    if(matcher.group(28) != null) {
                        String sequenceStr = matcher.group(28);
                        prediction.setSpectrum(spectra);
                        Sequence sequence = new Sequence();
                        sequence.setSequenceElement(new LinkedList<>());
                        sequence.setPeptideSequence(sequenceStr);
                        for (int i=0; i< sequenceStr.length();i++){
                            SEAminoAcid seAminoAcid = new SEAminoAcid();
                            AminoAcidFactory aminoAcidFactory = new AminoAcidFactory();
                            AminoAcid aminoAcid = aminoAcidFactory.createAminoAcid(String.valueOf(sequenceStr.charAt(i)));
                            seAminoAcid.setPositionInSequence(i + 2);
                            seAminoAcid.setAminoAcid(aminoAcid);
                            sequence.getSequenceElement().add(seAminoAcid);
                        }
                        prediction.setSequence(sequence);
                        targetDNMSO.getPrediction().add(prediction);
                    }
                    if(matcher.group(30) != null) {
                        Double highMSMSScore = Double.parseDouble(matcher.group(30));
                        scoreHighMSMSScore = new Score();
                        scoreHighMSMSScore.setScoreName("highMSMSScore");
                        scoreHighMSMSScore.setScoreValue(highMSMSScore);
                        prediction.getScore().add(scoreHighMSMSScore);
                    }
                    if(matcher.group(31) != null) {
                        Double highDiscriminantScore = Double.parseDouble(matcher.group(31));
                        scoreHighDiscriminantScore = new Score();
                        scoreHighDiscriminantScore.setScoreName("highDiscriminantScore");
                        scoreHighDiscriminantScore.setScoreValue(highDiscriminantScore);
                        prediction.getScore().add(scoreHighDiscriminantScore);
                    }
                    if(matcher.group(33) != null) {
                        Double pmtQualityScore = Double.parseDouble(matcher.group(33));
                        scorePMTQualityScore = new Score();
                        scorePMTQualityScore.setScoreName("pmtQualityScore");
                        scorePMTQualityScore.setScoreValue(pmtQualityScore);
                        prediction.getScore().add(scorePMTQualityScore);
                    }
                }
            }
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return targetDNMSO;
    }

    private boolean check(int peredictionSize, String size) {
        if(size!=null){
            return peredictionSize < Integer.parseInt(size);
        }
        return true;
    }

    @Override
    public boolean isValid(File file) {

        String line;
        try {
            BufferedReader in = new BufferedReader(new FileReader(file));
            while ((line = in.readLine()) != null) {
                if (line.startsWith("Sample Name")){
                    in.close();
                    return true;
                }
            }
            in.close();
        } catch (IOException e) {
            return false;
        }
        return false;
    }
}