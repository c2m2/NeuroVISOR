﻿using C2M2.NeuronalDynamics.Simulation;
using System;
using TMPro;
using UnityEngine;

[RequireComponent(typeof(SparseSolverTestv1))]
public class SimulationTimerLabel : MonoBehaviour
{
    public SparseSolverTestv1 sim = null;
    public TextMeshProUGUI timerText;

    /// <summary>
    /// Current time in simulation
    /// </summary>
    public double time;

    private void Start()
    {
        if (timerText == null) throw new LabelNotFoundException();
        timerText.text = time.ToString();
    }


    private void Update()
    {
        time = sim.GetSimulationTime();
        timerText.text = ToString();
    }
    public override string ToString()
    {
        if (time > 1000) return String.Format("{0:f0} s     {1:f0} ms", time/1000, time%1000);
        else return String.Format("{0:f0} ms", time);
    }

    public class LabelNotFoundException : Exception
    {
        public LabelNotFoundException() { }
        public LabelNotFoundException(string message) : base(message) { }
        public LabelNotFoundException(string message, Exception inner) : base(message, inner) { }
    }
}