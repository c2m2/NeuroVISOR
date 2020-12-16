﻿using C2M2.Simulation;
using System;
using System.Collections.Generic;
using TMPro;
using UnityEngine;

[RequireComponent(typeof(List<Canvas>))]
[RequireComponent(typeof(List<float>))]
public class RulerMeasure : MonoBehaviour
{
    public MeshSimulation sim = null;
    public List<Canvas> measurementDisplays;
    public List<float> numbers;
    private List<Tuple<TextMeshProUGUI, float>> markers = new List<Tuple<TextMeshProUGUI, float>>();
    private float relativeLength;
    private float initialRulerLength;
    private float minimumMarker = 0.05f; //earliest location on the ruler a marker can appear
    private float maximumMarker = 0.95f; //last location on the ruler a marker can appear
    private float markerSpacing = 0.05f; //minimum spacing between each marker

    // Start is called before the first frame update
    void Start()
    {
        relativeLength = 0;
        initialRulerLength = transform.lossyScale.z;
        CreateMarkers();
    }

    // Update is called once per frame
    void Update()
    {
        if (sim != null)
        {
            relativeLength = maximumMarker * (transform.lossyScale.z / sim.transform.localScale.z);

            int magnitude = GetMagnitude(relativeLength);
            string unit = GetUnit(magnitude);

            int siPrefixGroup = (int)Math.Floor(magnitude / 3.0);
            // length is a scaled version of relativelength so it is between 1 and 1000
            float length = (float)(relativeLength / Math.Pow(10, siPrefixGroup));
            UpdateMarkers(length, unit);
        }
    }

    private int GetMagnitude(float length)
    {
        double lengthLog10 = Math.Log10(length);
        return Convert.ToInt32(Math.Floor(lengthLog10));
    }

    private string GetUnit(int magnitude)
    {
        if (magnitude < -3)
        {
            // takes the magnitude and puts it in terms of nm by adding three. Then divides by 3 to get unit group and rounds down.
            int eTerm = (magnitude + 3) / 3;
            return " e" + 3 * eTerm + " nm";
        }
        else if (magnitude < 0) return " nm";
        else if (magnitude < 3) return " μm";
        else if (magnitude < 6) return " mm";
        else if (magnitude < 9) return " m";
        else if (magnitude <= 12) return " km";
        else if (magnitude > 12)
        {
            // takes the magnitude and puts it in terms of km by subtracting twelve. Then divides by 3 to get unit group and rounds down.
            int eTerm = (magnitude - 12) / 3;
            return " e" + 3 * eTerm + " km";
        }
        else
        {
            Debug.LogError("Invalid length inputted");
            return null;
        }
    }

    private void CreateMarkers()
    {
        foreach(Canvas measurementDisplay in measurementDisplays)
        {
            foreach(float number in numbers)
            {
                GameObject marker = new GameObject();
                marker.transform.SetParent(measurementDisplay.transform);
                TextMeshProUGUI markerText = marker.AddComponent<TextMeshProUGUI>();
                markerText.alignment = TextAlignmentOptions.Center;
                markerText.rectTransform.localPosition = new Vector3(0, 0, 0);
                markerText.fontSize = 0.03f;
                markerText.color = Color.black;
                marker.transform.localRotation = Quaternion.Euler(0,0,90);
                marker.transform.localScale = Vector3.one;
                markers.Add(new Tuple<TextMeshProUGUI, float>(markerText, number));
            }
            
        }
    }

    private void UpdateMarkers(float scaledRulerLength, string unit)
    {
        float previousLengthRatio = 0;
        foreach (Tuple<TextMeshProUGUI, float> marker in markers)
        {
            float markerNumber = marker.Item2;
            TextMeshProUGUI markerText = marker.Item1;
            float lengthRatio = (1/maximumMarker) * (markerNumber / scaledRulerLength);
            if (lengthRatio >= minimumMarker && lengthRatio <= maximumMarker && Math.Abs(lengthRatio-previousLengthRatio) > markerSpacing)
            {
                float rulerPoint = transform.localScale.z * (lengthRatio - .5f); // converts lengthRatio which goes from 0 to 1 to a point on the ruler which goes from -(length/2) to +(length/2)
                markerText.rectTransform.localPosition = new Vector3(rulerPoint, 0, 0);
                markerText.text = "― " + markerNumber + " " + unit + " ―";
                markerText.transform.localScale = new Vector3(markerText.transform.localScale.x, markerText.transform.localScale.y, transform.lossyScale.z/initialRulerLength);
                markerText.gameObject.SetActive(true);
                previousLengthRatio = lengthRatio;
            }
            else
            {
                markerText.gameObject.SetActive(false);
            }
        }
    }

}
