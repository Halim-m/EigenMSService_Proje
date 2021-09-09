//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vhudson-jaxb-ri-2.1-833 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2012.08.07 at 10:34:54 AM EEST 
//


package domain.pepxml;

import javax.xml.bind.annotation.*;


/**
 * <p>Java class for anonymous complex type. <p>The following schema fragment specifies the expected content contained within this class. <pre> &lt;complexType> &lt;complexContent> &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType"> &lt;sequence> &lt;element ref="{http://regis-web.systemsbiology.net/pepXML}asapratio_peptide_data"/> &lt;/sequence> &lt;attribute name="mean" use="required" type="{http://www.w3.org/2001/XMLSchema}float" /> &lt;attribute name="error" use="required" type="{http://www.w3.org/2001/XMLSchema}float" /> &lt;attribute name="heavy2light_mean" use="required" type="{http://www.w3.org/2001/XMLSchema}float" /> &lt;attribute name="heavy2light_error" use="required" type="{http://www.w3.org/2001/XMLSchema}float" /> &lt;/restriction> &lt;/complexContent> &lt;/complexType> </pre>
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "asapratioPeptideData"
})
@XmlRootElement(name = "asapratio_result")
public class AsapratioResult {

    /**
	 * @uml.property  name="asapratioPeptideData"
	 * @uml.associationEnd  
	 */
    @XmlElement(name = "asapratio_peptide_data", required = true)
    protected AsapratioPeptideData asapratioPeptideData;
    /**
	 * @uml.property  name="mean"
	 */
    @XmlAttribute(required = true)
    protected float mean;
    /**
	 * @uml.property  name="error"
	 */
    @XmlAttribute(required = true)
    protected float error;
    /**
	 * @uml.property  name="heavy2LightMean"
	 */
    @XmlAttribute(name = "heavy2light_mean", required = true)
    protected float heavy2LightMean;
    /**
	 * @uml.property  name="heavy2LightError"
	 */
    @XmlAttribute(name = "heavy2light_error", required = true)
    protected float heavy2LightError;

    /**
	 * Gets the value of the asapratioPeptideData property.
	 * @return  possible object is {@link AsapratioPeptideData  }
	 * @uml.property  name="asapratioPeptideData"
	 */
    public AsapratioPeptideData getAsapratioPeptideData() {
        return asapratioPeptideData;
    }

    /**
	 * Sets the value of the asapratioPeptideData property.
	 * @param value  allowed object is {@link AsapratioPeptideData  }
	 * @uml.property  name="asapratioPeptideData"
	 */
    public void setAsapratioPeptideData(AsapratioPeptideData value) {
        this.asapratioPeptideData = value;
    }

    /**
	 * Gets the value of the mean property.
	 * @uml.property  name="mean"
	 */
    public float getMean() {
        return mean;
    }

    /**
	 * Sets the value of the mean property.
	 * @uml.property  name="mean"
	 */
    public void setMean(float value) {
        this.mean = value;
    }

    /**
	 * Gets the value of the error property.
	 * @uml.property  name="error"
	 */
    public float getError() {
        return error;
    }

    /**
	 * Sets the value of the error property.
	 * @uml.property  name="error"
	 */
    public void setError(float value) {
        this.error = value;
    }

    /**
	 * Gets the value of the heavy2LightMean property.
	 * @uml.property  name="heavy2LightMean"
	 */
    public float getHeavy2LightMean() {
        return heavy2LightMean;
    }

    /**
	 * Sets the value of the heavy2LightMean property.
	 * @uml.property  name="heavy2LightMean"
	 */
    public void setHeavy2LightMean(float value) {
        this.heavy2LightMean = value;
    }

    /**
	 * Gets the value of the heavy2LightError property.
	 * @uml.property  name="heavy2LightError"
	 */
    public float getHeavy2LightError() {
        return heavy2LightError;
    }

    /**
	 * Sets the value of the heavy2LightError property.
	 * @uml.property  name="heavy2LightError"
	 */
    public void setHeavy2LightError(float value) {
        this.heavy2LightError = value;
    }

}
