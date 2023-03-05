#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "motion.h"
#include "interpolator.h"
#include "types.h"
#include "transform.h"
#include "performanceCounter.h"
#include <stdexcept>
#include <vector>

// Initialize performance counter to measure interpolation time
PerformanceCounter pc;

/* Generate random intervals between keyframes.
 * N: The maximum possible interval length.
 * totalLength: The length of input video in frames.
 * return: A std::vector that represents **intervals** between keyframes.
 * 
 * Note: The first and last frame are always keyframes.
 * */
std::vector<int> randomKeyframes(int N, int totalLength)
{
  srand(0);  // set random seed
  std::vector<int> keyframeIntervals;  // vector to hold the intervals
  totalLength -= 1;
  int currentFrame = 0;
  while (currentFrame < totalLength)
  {
    int maxInterval = fmin(N, totalLength - currentFrame);
    int interval = rand() % maxInterval + 1;
    keyframeIntervals.push_back(interval);
    currentFrame += interval + 1;
  }
  int excess = currentFrame - totalLength;
  if (excess > 0)
  {
    keyframeIntervals[keyframeIntervals.size() - 1] -= excess;
  }
  return keyframeIntervals;
}

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

  // Perform the interpolation on uniform keyframe intervals
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER) && N >= 0)
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION) && N >= 0)
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER) && N >= 0)
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION) && N >= 0)
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  // Perform interpolation on non-uniform keyframe intervals (randomly generated)
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER) && N < 0)
    LinearInterpolationEulerRandom(pInputMotion, *pOutputMotion, -N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION) && N < 0)
    LinearInterpolationQuaternionRandom(pInputMotion, *pOutputMotion, -N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER) && N < 0)
    BezierInterpolationEulerRandom(pInputMotion, *pOutputMotion, -N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION) && N < 0)
    BezierInterpolationQuaternionRandom(pInputMotion, *pOutputMotion, -N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  // Starts timing for performance analysis
  pc.StartCounter();
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::LinearInterpolationEulerRandom(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  // Generates random keyframe intervals
  std::vector<int> keyframeIntervals = randomKeyframes(N, inputLength);
  // Prints the random intervals
  printf("The generated random keyframe intervals are: [");
  for (int i = 0; i < keyframeIntervals.size(); ++i)
  {
    printf("%d", keyframeIntervals[i]);
    if (i < keyframeIntervals.size() - 1)
      printf(", ");
  }
  printf("]\n");

  int startKeyframe = 0;
  // Starts timing for performance analysis
  pc.StartCounter();
  for (int N: keyframeIntervals)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON) 
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  } 
  else 
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}

void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
  double Rx[4][4], Ry[4][4], Rz[4][4], temp[4][4], result[4][4];
  /* Creates rotation matrices along x, y, z axes, respectively */
  rotationX(Rx, angles[0]);
  rotationY(Ry, angles[1]);
  rotationZ(Rz, angles[2]);
  /* Forms the rotation matrix */
  matrix_mult(Rz, Ry, temp);  // Ry = Rz * Ry
  matrix_mult(temp, Rx, result);  // result = Rz * Ry * Rx
  /* Copies the value from (homogeneous) rotation matrix (4x4) to output rotation matrix (3x3) */
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      R[j + 3 * i] = result[i][j];
}

void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  if (N == inputLength - 2)
  {
    printf("Reducing to `LinearInterpolationEuler` since `N=%d` while `inputLength=%d`.\n", N, inputLength);
    return LinearInterpolationEuler(pInputMotion, pOutputMotion, N);
  }

  int startKeyframe = 0;
  // Starts timing for performance analysis
  pc.StartCounter();
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;
    int prevKeyframe = startKeyframe - N - 1;
    int nextKeyframe = endKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
    Posture * prevPosture;
    Posture * nextPosture;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    /* Initializes for Bezier Interpolation */
    vector posPrev, posStart, posEnd, posNext, temp;
    posStart = startPosture->root_pos;
    posEnd = endPosture->root_pos;
    if (prevKeyframe >= 0)
    {
      prevPosture = pInputMotion->GetPosture(prevKeyframe);
      posPrev = prevPosture->root_pos;
    }
    if (nextKeyframe < inputLength)
    {
      nextPosture = pInputMotion->GetPosture(nextKeyframe);
      posNext = nextPosture->root_pos;
    }

    vector posA, posB;  // control points for root position Bezier interpolation
    if (startKeyframe == 0)  // The first key frame
    {
      // [Special case] Calculates posA
      temp = lerp(2.0, posNext, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    else if (nextKeyframe > inputLength - 1)  // The second last keyframe, i.e., no keyframes after current end keyframe
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // [Special case] Calculates posB
      temp = lerp(2.0, posPrev, posStart);
      posB = lerp(1.0 / 3, posEnd, temp);
    }
    else  // In-between keyframes
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    vector rotStart[MAX_BONES_IN_ASF_FILE];
    vector rotEnd[MAX_BONES_IN_ASF_FILE];
    vector rotA[MAX_BONES_IN_ASF_FILE];  // control points A for bones
    vector rotB[MAX_BONES_IN_ASF_FILE];  // control points B for bones
    vector rotPrev, rotNext, rotTemp;
    for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
    {
      rotStart[bone] = startPosture->bone_rotation[bone];
      rotEnd[bone] = endPosture->bone_rotation[bone];
      if (startKeyframe == 0)
      {
        rotNext = nextPosture->bone_rotation[bone];
        // [Special case] Calculates rotA for each bone
        rotTemp = lerp(2.0, rotNext, rotEnd[bone]);
        rotA[bone] = lerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = lerp(2.0, rotStart[bone], rotEnd[bone]);
        rotTemp = lerp(0.5, rotTemp, rotNext);
        rotB[bone] = lerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
      else if (nextKeyframe > inputLength - 1)
      {
        rotPrev = prevPosture->bone_rotation[bone];
        // Calculates rotA for each bone
        rotTemp = lerp(2.0, rotPrev, rotStart[bone]);
        rotTemp = lerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = lerp(1.0 / 3, rotStart[bone], rotTemp);
        // [Special case] Calculates rotB for each bone
        rotTemp = lerp(2.0, rotPrev, rotStart[bone]);
        rotB[bone] = lerp(1.0 / 3, rotEnd[bone], rotTemp);
      }
      else  // In-between keyframes
      {
        rotPrev = prevPosture->bone_rotation[bone];
        rotNext = nextPosture->bone_rotation[bone];
        // Calculates rotA for each bone
        rotTemp = lerp(2.0, rotPrev, rotStart[bone]);
        rotTemp = lerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = lerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = lerp(2.0, rotStart[bone], rotEnd[bone]);
        rotTemp = lerp(0.5, rotTemp, rotNext);
        rotB[bone] = lerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
    }

    // interpolate in between
    for (int frame = 1; frame <= N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N + 1);

      // interpolate root position in Bezier Euler
      interpolatedPosture.root_pos = DeCasteljauEuler(t, posStart, posA, posB, posEnd);

      // interpolate bone rotations in Bezier Euler
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, rotStart[bone], rotA[bone], rotB[bone], rotEnd[bone]);

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for (int frame = startKeyframe + 1; frame < inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::BezierInterpolationEulerRandom(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  // Generates random keyframe intervals
  std::vector<int> keyframeIntervals = randomKeyframes(N, inputLength);
  // Prints the random intervals
  printf("The generated random keyframe intervals are: [");
  for (int i = 0; i < keyframeIntervals.size(); ++i)
  {
    printf("%d", keyframeIntervals[i]);
    if (i < keyframeIntervals.size() - 1)
      printf(", ");
  }
  printf("]\n");

  int startKeyframe = 0;
  // Starts timing for performance analysis
  pc.StartCounter();
  for (int i = 0; i < keyframeIntervals.size(); ++i)
  {
    N = keyframeIntervals[i];
    int endKeyframe = startKeyframe + N + 1;
    int prevKeyframe = -1;
    int nextKeyframe = inputLength;
    try {
      prevKeyframe = startKeyframe - keyframeIntervals.at(i - 1) - 1;
    }
    catch (const std::out_of_range& e) { }
    try {
      nextKeyframe = endKeyframe + keyframeIntervals.at(i + 1) + 1;
    } catch (const std::out_of_range& e) { }

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
    Posture * prevPosture;
    Posture * nextPosture;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    /* Initializes for Bezier Interpolation */
    vector posPrev, posStart, posEnd, posNext, temp;
    posStart = startPosture->root_pos;
    posEnd = endPosture->root_pos;
    if (prevKeyframe >= 0)
    {
      prevPosture = pInputMotion->GetPosture(prevKeyframe);
      posPrev = prevPosture->root_pos;
    }
    if (nextKeyframe < inputLength)
    {
      nextPosture = pInputMotion->GetPosture(nextKeyframe);
      posNext = nextPosture->root_pos;
    }

    vector posA, posB;  // control points for root position Bezier interpolation
    if (startKeyframe == 0)  // The first key frame
    {
      // [Special case] Calculates posA
      temp = lerp(2.0, posNext, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    else if (nextKeyframe > inputLength - 1)  // The second last keyframe, i.e., no keyframes after current end keyframe
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // [Special case] Calculates posB
      temp = lerp(2.0, posPrev, posStart);
      posB = lerp(1.0 / 3, posEnd, temp);
    }
    else  // In-between keyframes
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    vector rotStart[MAX_BONES_IN_ASF_FILE];
    vector rotEnd[MAX_BONES_IN_ASF_FILE];
    vector rotA[MAX_BONES_IN_ASF_FILE];  // control points A for bones
    vector rotB[MAX_BONES_IN_ASF_FILE];  // control points B for bones
    vector rotPrev, rotNext, rotTemp;
    for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
    {
      rotStart[bone] = startPosture->bone_rotation[bone];
      rotEnd[bone] = endPosture->bone_rotation[bone];
      if (startKeyframe == 0)
      {
        rotNext = nextPosture->bone_rotation[bone];
        // [Special case] Calculates rotA for each bone
        rotTemp = lerp(2.0, rotNext, rotEnd[bone]);
        rotA[bone] = lerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = lerp(2.0, rotStart[bone], rotEnd[bone]);
        rotTemp = lerp(0.5, rotTemp, rotNext);
        rotB[bone] = lerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
      else if (nextKeyframe > inputLength - 1)
      {
        rotPrev = prevPosture->bone_rotation[bone];
        // Calculates rotA for each bone
        rotTemp = lerp(2.0, rotPrev, rotStart[bone]);
        rotTemp = lerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = lerp(1.0 / 3, rotStart[bone], rotTemp);
        // [Special case] Calculates rotB for each bone
        rotTemp = lerp(2.0, rotPrev, rotStart[bone]);
        rotB[bone] = lerp(1.0 / 3, rotEnd[bone], rotTemp);
      }
      else  // In-between keyframes
      {
        rotPrev = prevPosture->bone_rotation[bone];
        rotNext = nextPosture->bone_rotation[bone];
        // Calculates rotA for each bone
        rotTemp = lerp(2.0, rotPrev, rotStart[bone]);
        rotTemp = lerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = lerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = lerp(2.0, rotStart[bone], rotEnd[bone]);
        rotTemp = lerp(0.5, rotTemp, rotNext);
        rotB[bone] = lerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
    }

    // interpolate in between
    for (int frame = 1; frame <= N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N + 1);

      // interpolate root position in Bezier Euler
      interpolatedPosture.root_pos = DeCasteljauEuler(t, posStart, posA, posB, posEnd);

      // interpolate bone rotations in Bezier Euler
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, rotStart[bone], rotA[bone], rotB[bone], rotEnd[bone]);

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for (int frame = startKeyframe + 1; frame < inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  Quaternion<double> qStart, qEnd, qResult;
  // Starts timing for performance analysis
  pc.StartCounter();
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations in SLERP
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
        Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd);
        qResult = Slerp(t, qStart, qEnd);
        Quaternion2Euler(qResult, interpolatedPosture.bone_rotation[bone].p);
      }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::LinearInterpolationQuaternionRandom(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  // Generates random keyframe intervals
  std::vector<int> keyframeIntervals = randomKeyframes(N, inputLength);
  // Prints the random intervals
  printf("The generated random keyframe intervals are: [");
  for (int i = 0; i < keyframeIntervals.size(); ++i)
  {
    printf("%d", keyframeIntervals[i]);
    if (i < keyframeIntervals.size() - 1)
      printf(", ");
  }
  printf("]\n");

  int startKeyframe = 0;
  Quaternion<double> qStart, qEnd, qResult;

  // Starts timing for performance analysis
  pc.StartCounter();
  for (int N: keyframeIntervals)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations in SLERP
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
        Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd);
        qResult = Slerp(t, qStart, qEnd);
        Quaternion2Euler(qResult, interpolatedPosture.bone_rotation[bone].p);
      }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  if (N == inputLength - 2)
  {
    printf("Reducing to `LinearInterpolationQuaternion` since `N=%d` while `inputLength=%d`.\n", N, inputLength);
    return LinearInterpolationQuaternion(pInputMotion, pOutputMotion, N);
  }

  int startKeyframe = 0;
  // Starts timing for performance analysis
  pc.StartCounter();
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;
    int prevKeyframe = startKeyframe - N - 1;
    int nextKeyframe = endKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
    Posture * prevPosture;
    Posture * nextPosture;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    /* Initializes for Bezier Interpolation */
    vector posPrev, posStart, posEnd, posNext, temp;
    posStart = startPosture->root_pos;
    posEnd = endPosture->root_pos;
    if (prevKeyframe >= 0)
    {
      prevPosture = pInputMotion->GetPosture(prevKeyframe);
      posPrev = prevPosture->root_pos;
    }
    if (nextKeyframe < inputLength)
    {
      nextPosture = pInputMotion->GetPosture(nextKeyframe);
      posNext = nextPosture->root_pos;
    }

    vector posA, posB;  // control points for root position Bezier interpolation
    if (startKeyframe == 0)  // The first key frame
    {
      // [Special case] Calculates posA
      temp = lerp(2.0, posNext, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    else if (nextKeyframe > inputLength - 1)  // The second last keyframe, i.e., no keyframes after current end keyframe
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // [Special case] Calculates posB
      temp = lerp(2.0, posPrev, posStart);
      posB = lerp(1.0 / 3, posEnd, temp);
    }
    else  // In-between keyframes
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    Quaternion<double> rotStart[MAX_BONES_IN_ASF_FILE];
    Quaternion<double> rotEnd[MAX_BONES_IN_ASF_FILE];
    Quaternion<double> rotA[MAX_BONES_IN_ASF_FILE];  // control points A for bones
    Quaternion<double> rotB[MAX_BONES_IN_ASF_FILE];  // control points B for bones
    Quaternion<double> rotPrev, rotNext, rotTemp;
    for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
    {
      Euler2Quaternion(startPosture->bone_rotation[bone].p, rotStart[bone]);
      Euler2Quaternion(endPosture->bone_rotation[bone].p, rotEnd[bone]);
      if (startKeyframe == 0)
      {
        Euler2Quaternion(nextPosture->bone_rotation[bone].p, rotNext);
        // [Special case] Calculates rotA for each bone
        rotTemp = Double(rotNext, rotEnd[bone]);
        rotA[bone] = Slerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = Double(rotStart[bone], rotEnd[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotNext);
        rotB[bone] = Slerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
      else if (nextKeyframe > inputLength - 1)
      {
        Euler2Quaternion(prevPosture->bone_rotation[bone].p, rotPrev);
        // Calculates rotA for each bone
        rotTemp = Double(rotPrev, rotStart[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = Slerp(1.0 / 3, rotStart[bone], rotTemp);
        // [Special case] Calculates rotB for each bone
        rotTemp = Double(rotPrev, rotStart[bone]);
        rotB[bone] = Slerp(1.0 / 3, rotEnd[bone], rotTemp);
      }
      else  // In-between keyframes
      {
        Euler2Quaternion(prevPosture->bone_rotation[bone].p, rotPrev);
        Euler2Quaternion(nextPosture->bone_rotation[bone].p, rotNext);
        // Calculates rotA for each bone
        rotTemp = Double(rotPrev, rotStart[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = Slerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = Double(rotStart[bone], rotEnd[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotNext);
        rotB[bone] = Slerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
    }

    // interpolate in between
    Quaternion<double> qResult;
    for (int frame = 1; frame <= N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N + 1);

      // interpolate root position in Bezier Euler
      interpolatedPosture.root_pos = DeCasteljauEuler(t, posStart, posA, posB, posEnd);

      // interpolate bone rotations in Bezier Quaternion
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        qResult = DeCasteljauQuaternion(t, rotStart[bone], rotA[bone], rotB[bone], rotEnd[bone]);
        Quaternion2Euler(qResult, interpolatedPosture.bone_rotation[bone].p);
      }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for (int frame = startKeyframe + 1; frame < inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::BezierInterpolationQuaternionRandom(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  // Generates random keyframe intervals
  std::vector<int> keyframeIntervals = randomKeyframes(N, inputLength);
  // Prints the random intervals
  printf("The generated random keyframe intervals are: [");
  for (int i = 0; i < keyframeIntervals.size(); ++i)
  {
    printf("%d", keyframeIntervals[i]);
    if (i < keyframeIntervals.size() - 1)
      printf(", ");
  }
  printf("]\n");

  int startKeyframe = 0;
  // Starts timing for performance analysis
  pc.StartCounter();
  for (int i = 0; i < keyframeIntervals.size(); ++i)
  {
    N = keyframeIntervals[i];
    int endKeyframe = startKeyframe + N + 1;
    int prevKeyframe = -1;
    int nextKeyframe = inputLength;
    try {
      prevKeyframe = startKeyframe - keyframeIntervals.at(i - 1) - 1;
    }
    catch (const std::out_of_range& e) { }
    try {
      nextKeyframe = endKeyframe + keyframeIntervals.at(i + 1) + 1;
    } catch (const std::out_of_range& e) { }

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
    Posture * prevPosture;
    Posture * nextPosture;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    /* Initializes for Bezier Interpolation */
    vector posPrev, posStart, posEnd, posNext, temp;
    posStart = startPosture->root_pos;
    posEnd = endPosture->root_pos;
    if (prevKeyframe >= 0)
    {
      prevPosture = pInputMotion->GetPosture(prevKeyframe);
      posPrev = prevPosture->root_pos;
    }
    if (nextKeyframe < inputLength)
    {
      nextPosture = pInputMotion->GetPosture(nextKeyframe);
      posNext = nextPosture->root_pos;
    }

    vector posA, posB;  // control points for root position Bezier interpolation
    if (startKeyframe == 0)  // The first key frame
    {
      // [Special case] Calculates posA
      temp = lerp(2.0, posNext, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    else if (nextKeyframe > inputLength - 1)  // The second last keyframe, i.e., no keyframes after current end keyframe
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // [Special case] Calculates posB
      temp = lerp(2.0, posPrev, posStart);
      posB = lerp(1.0 / 3, posEnd, temp);
    }
    else  // In-between keyframes
    {
      // Calculates posA
      temp = lerp(2.0, posPrev, posStart);
      temp = lerp(0.5, temp, posEnd);
      posA = lerp(1.0 / 3, posStart, temp);
      // Calculates posB
      temp = lerp(2.0, posStart, posEnd);
      temp = lerp(0.5, temp, posNext);
      posB = lerp(-1.0 / 3, posEnd, temp);
    }
    Quaternion<double> rotStart[MAX_BONES_IN_ASF_FILE];
    Quaternion<double> rotEnd[MAX_BONES_IN_ASF_FILE];
    Quaternion<double> rotA[MAX_BONES_IN_ASF_FILE];  // control points A for bones
    Quaternion<double> rotB[MAX_BONES_IN_ASF_FILE];  // control points B for bones
    Quaternion<double> rotPrev, rotNext, rotTemp;
    for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
    {
      Euler2Quaternion(startPosture->bone_rotation[bone].p, rotStart[bone]);
      Euler2Quaternion(endPosture->bone_rotation[bone].p, rotEnd[bone]);
      if (startKeyframe == 0)
      {
        Euler2Quaternion(nextPosture->bone_rotation[bone].p, rotNext);
        // [Special case] Calculates rotA for each bone
        rotTemp = Double(rotNext, rotEnd[bone]);
        rotA[bone] = Slerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = Double(rotStart[bone], rotEnd[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotNext);
        rotB[bone] = Slerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
      else if (nextKeyframe > inputLength - 1)
      {
        Euler2Quaternion(prevPosture->bone_rotation[bone].p, rotPrev);
        // Calculates rotA for each bone
        rotTemp = Double(rotPrev, rotStart[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = Slerp(1.0 / 3, rotStart[bone], rotTemp);
        // [Special case] Calculates rotB for each bone
        rotTemp = Double(rotPrev, rotStart[bone]);
        rotB[bone] = Slerp(1.0 / 3, rotEnd[bone], rotTemp);
      }
      else  // In-between keyframes
      {
        Euler2Quaternion(prevPosture->bone_rotation[bone].p, rotPrev);
        Euler2Quaternion(nextPosture->bone_rotation[bone].p, rotNext);
        // Calculates rotA for each bone
        rotTemp = Double(rotPrev, rotStart[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotEnd[bone]);
        rotA[bone] = Slerp(1.0 / 3, rotStart[bone], rotTemp);
        // Calculates rotB for each bone
        rotTemp = Double(rotStart[bone], rotEnd[bone]);
        rotTemp = Slerp(0.5, rotTemp, rotNext);
        rotB[bone] = Slerp(-1.0 / 3, rotEnd[bone], rotTemp);
      }
    }

    // interpolate in between
    Quaternion<double> qResult;
    for (int frame = 1; frame <= N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N + 1);

      // interpolate root position in Bezier Euler
      interpolatedPosture.root_pos = DeCasteljauEuler(t, posStart, posA, posB, posEnd);

      // interpolate bone rotations in Bezier Quaternion
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        qResult = DeCasteljauQuaternion(t, rotStart[bone], rotA[bone], rotB[bone], rotEnd[bone]);
        Quaternion2Euler(qResult, interpolatedPosture.bone_rotation[bone].p);
      }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  // Stop timing for performance analysis
  pc.StopCounter();
  printf("Interpolation time: %f ms.\n", pc.GetElapsedTime());
  for (int frame = startKeyframe + 1; frame < inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) 
{
  double R[9];
  Euler2Rotation(angles, R);
  q = Quaternion<double>::Matrix2Quaternion(R);
  q.Normalize();
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) 
{
  double R[9];
  q.Quaternion2Matrix(R);
  Rotation2Euler(R, angles);
}

/* Linear interpolation on 3-dim vectors */
vector lerp(double t, vector & start, vector & end)
{
  vector c = (end - start) * t + start;
  return c;
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd_)
{
  Quaternion<double> result;
  double cosTheta, theta;
  cosTheta = qStart.Gets() * qEnd_.Gets() + qStart.Getx() * qEnd_.Getx() + qStart.Gety() * qEnd_.Gety() + qStart.Getz() * qEnd_.Getz();
  cosTheta = clamp(cosTheta, -1.0, 1.0);  // avoid NaN problem
  if (cosTheta < 0)  // flip to the shortest arc
  {
    qEnd_ = -1.0 * qEnd_;
    cosTheta *= -1.0;
  }

  theta = acos(cosTheta);

  /* Edge case: qStart and qEnd_ represent the same rotation */
  if (theta == 0.0)
    return qStart;

  /* General case */
  result = sin((1.0 - t) * theta) / sin(theta) * qStart + sin(t * theta) / sin(theta) * qEnd_;
  return result;
}

Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
  Quaternion<double> result;
  result = Slerp(2.0, p, q);
  return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
  vector result, Q0, Q1, Q2, R0, R1;
  Q0 = p0 * (1 - t) + p1 * t;
  Q1 = p1 * (1 - t) + p2 * t;
  Q2 = p2 * (1 - t) + p3 * t;
  R0 = Q0 * (1 - t) + Q1 * t;
  R1 = Q1 * (1 - t) + Q2 * t;
  result = R0 * (1 - t) + R1 * t;
  return result;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
  Quaternion<double> result, Q0, Q1, Q2, R0, R1;
  Q0 = Slerp(t, p0, p1);
  Q1 = Slerp(t, p1, p2);
  Q2 = Slerp(t, p2, p3);
  R0 = Slerp(t, Q0, Q1);
  R1 = Slerp(t, Q1, Q2);
  result = Slerp(t, R0, R1);
  return result;
}

