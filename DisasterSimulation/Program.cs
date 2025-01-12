using System;
using System.Text.Json;
using DisasterSimulation;

DirectoryInfo directoryInfo = new("../../../../source/faces");
IEnumerable<FileInfo> files = directoryInfo.EnumerateFiles("*", SearchOption.AllDirectories);
foreach (FileInfo file in files)
{
    StreamReader streamReader = new(file.FullName);
    string json = streamReader.ReadToEnd();
    streamReader.Close();
    FaceData[]? faceDataPart = JsonSerializer.Deserialize<FaceData[]>(json);
    Console.WriteLine(json);
    //Console.WriteLine(faceDataPart?[0].triangle[0].x.ToString());
    //Console.WriteLine(faceDataPart?[0].index.ToString());
}

DisasterSimulation.Simulator simulator = new([]);
simulator.Start(0);