﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using C2M2.Interaction.UI;
using C2M2.Utils;
namespace C2M2.Interaction.VR
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// <summary> Add movement functionality to Oculus HMD headset emulator,
    ///           add some additional small functionality like disabling hands in emulator mode,
    ///           allow slow movement </summary>
    ///
    /// <remarks>   Jacob Wells, 4/30/2020. </remarks>
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    [RequireComponent(typeof(MovementController))]
    public class MovingOVRHeadsetEmulator : OVRHeadsetEmulator
    {
        public KeyCode[] slowMoveKeys = new KeyCode[] { KeyCode.LeftShift, KeyCode.RightShift };

        public float speed = 0.1f;
        public float slowSpeed = 0.025f;
        private MovementController controls = null;
        private bool SlowMoving
        {
            get
            {
                bool slowMoving = false;
                foreach (KeyCode key in slowMoveKeys)
                {
                    if (Input.GetKey(key))
                        slowMoving = true;
                }
                return slowMoving;
            }
        }

        private void Awake()
        {
            controls = GetComponent<MovementController>() ?? gameObject.AddComponent<MovementController>();
            controls.speed = speed;
        }

        private void OnEnable()
        {
            StartCoroutine(CheckSpeed());

            if (controls != null) controls.enabled = true;
        }

        private void OnDisable()
        {
            StopCoroutine(CheckSpeed());

            if (controls != null) controls.enabled = false;
        }

        private IEnumerator CheckSpeed()
        {
            while (true)
            {
                controls.speed = SlowMoving ? slowSpeed : speed;
                yield return null;
            }
        }

       
    }
}
